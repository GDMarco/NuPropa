#include "NeutrinoNeutrinoInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <string>

namespace crpropa {

// The parent's constructor need to be called on initialization!
NeutrinoNeutrinoInteraction::NeutrinoNeutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, bool haveSecondaries, double limit) : Module() {
    setNeutrinoField(neutrinoField);
    setHaveSecondary(haveSecondaries);
    setLimit(limit);
}

void NeutrinoNeutrinoInteraction::setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField) {
    this->neutrinoField = neutrinoField;
    this->neutrinoFieldID = neutrinoField->getParticleID(); 
    std::string fname = neutrinoField->getFieldName();
    setDescription("NeutrinoNeutrinoInteraction::Module" + fname);
    
    fileNuNu = getDataPath("NeutrinoNeutrinoInteraction/NeutrinoNeutrinoElastic/rate_" + fname + ".txt");
    fileNuiNuj = getDataPath("NeutrinoNeutrinoInteraction/NeutrinoiNeutrinojElastic/rate_" + fname + ".txt");
    fileNuiNuj = getDataPath();
    initRate(fileNuNu, fileNuiNuj);
}

void NeutrinoNeutrinoInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoNeutrinoInteraction::setLimit(double limit) {
    this->limit = limit;
}


void NeutrinoNeutrinoInteraction::initRate(std::string fileNuNu, std::string fileNuiNuj) {
    std::ifstream infileNuNu(fileNuNu.c_str());
    std::ifstream infileNuiNuj(fileNuiNuj.c_str());
    
    tabEnergy.clear();
    tabRate.clear();
    
    if (!infileNuNu.good())
        throw std::runtime_error("NeutrinoNeutrinoInteraction: could not open file" + filename);
    
    std::vector<double> vecEnergyNuNu;
    std::vector<double> vecRateNuNu;
    
    while (infileNuNu.good()) {
        if (infileNuNu.peek() != '#') {
            double a, b;
            infile >> a >> b;
            if (infileNuNu) {
                vecEnergyNuNu.push_back(pow(10, a) * eV);
                vecRateNuNu.push_back(b / Mpc);
            }
        }
        infileNuNu.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infileNuNu.close();
    
    this->tabEnergy.push_back(vecEnergyNuNu);
    this->tabRate.push_back(vecRateNuNu);
    
    if (!infileNuiNuj.good())
        throw std::runtime_error("NeutrinoiNeutrinojInteraction: could not open file" + filename);
    
    std::vector<double> vecEnergyNuiNuj;
    std::vector<double> vecRateNuiNuj;
    
    while (infileNuiNuj.good()) {
        if (infileNuiNuj.peek() != '#') {
            double a, b;
            infile >> a >> b;
            if (infileNuiNuj) {
                vecEnergyNuiNuj.push_back(pow(10, a) * eV);
                vecRateNuiNuj.push_back(b / Mpc);
            }
        }
        infileNuiNuj.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infileNuiNuj.close();
    
    this->tabEnergy.push_back(vecEnergyNuiNuj);
    this->tabRate.push_back(vecRateNuiNuj);
}


void NeutrinoNeutrinoInteraction::performInteraction(Candidate *candidate) const {
    
    candidate->setActive(false);
    
    if (not haveSecondaries)
        return;
    
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    double ID = candidate->current.getID();
    double w = 1; // no thinning, neither useful
    
    double Enu = E / 10.; // naive value, waiting for the differential cross sections
    double Enu2 = E - Enu;
    Random &random = Random::instance();
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    
    if (haveSecondaries)
        candidate->addSecondary(ID, Enu / (1 + z), pos, w, interactionTag);
        candidate->addSecondary(this->neutrinoFieldID, Enu2 / (1 +z), pos, w, interactionTag);
}


void NeutrinoNeutrinoInteraction::process(Candidate *candidate) const
{
    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = (1 + z) * candidate->current.getEnergy();
    double ID = candidate->current.getID();
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;
    
    std::vector<double> vecEnergy;
    std::vector<double> vecRate;
    
    if (ID == this->neutrinoFieldID) {
        vecEnergy = this->tabEnergy[0];
        vecRate = this->tabRate[0];
        setInteractionTag("NuNuInt");
    } else {
        vecEnergy = this->tabEnergy[1];
        vecRate = this->tabRate[1];
        setInteractionTag("NuiNujInt");
    }
    
    // check if in tabulated energy range
    if (E < vecEnergy.front() or (E > vecEnergy.back()))
        return;
        
    // interaction rate
    double rate = interpolate(E, vecEnergy, vecRate);
    rate *= pow_integer<2>(1 + z) * neutrinoField->getRedshiftScaling(z); // to check how it defined this function in NeutrinoBackground.h
        
    // check for interaction
    Random &random = Random::instance();
    double randDistance = -log(random.rand()) / rate;
    double step = candidate->getCurrentStep();
    if (step < randDistance) {
        candidate->limitNextStep(limit / rate);
        return;
    } else { // after performing interaction neutrino ceases to exist (hence return)
        performInteraction(candidate);
        return;
    }
}

void NeutrinoNeutrinoInteraction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string NeutrinoPhotonInteraction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace crpropa

