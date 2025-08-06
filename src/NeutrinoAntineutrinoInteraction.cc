#include "nupropa/NeutrinoAntineutrinoInteraction.h"
#include "nupropa/Channels.h"
#include "nupropa/ChannelsBundle.h"
#include "nupropa/NeutrinoField.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/ParticleData.h"
#include "nupropa/NeutrinoMixing.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include <string>
#include <filesystem>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <cmath>

namespace nupropa {

using namespace crpropa;

NeutrinoAntineutrinoInteraction::NeutrinoAntineutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, ref_ptr<Channels> channels, bool haveSecondaries, double limit, ref_ptr<NeutrinoMixing> neutrinoMixing) : Module() { // double thinning
    
    setChannels(channels);
    setNeutrinoField(neutrinoField);
    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    //setThinning(thinning);
    setNeutrinoMixing(neutrinoMixing);
    
}

/**
 One assumes that the tables have been produced within the same computation, so the energy bins are the same.
 
 In the case of the neutrino-antineutrino interactions the elastic scatterings furnish the tabEnergy, since they do not have interaction energy thresholds (almost)!
 */

// channels have an unordered map <bool, "interaction">, maybe interaction tag? setChannels should be a way of filling this object, knowing how they are ordered!
void NeutrinoAntineutrinoInteraction::setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField) {
    this->neutrinoField = neutrinoField;
    this->neutrinoFieldMass = neutrinoField->getMass();
    std::string fname = neutrinoField->getFieldName();
    setDescription("NeutrinoAntineutrinoInteraction::Module" + fname);
    setInteractionTag("NuNuInt");
    
    setChannelsBundle(this->channels, fname);
    // initRate(fname);
}

void NeutrinoAntineutrinoInteraction::setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    this->neutrinoMixing = neutrinoMixing;
}

void NeutrinoAntineutrinoInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoAntineutrinoInteraction::setLimit(double limit) {
    this->limit = limit;
}

/**
void NeutrinoNeutrinoInteraction::setThinning(double thinning) {
    this->thinning = thinning;
}
*/

 void NeutrinoAntineutrinoInteraction::setChannels(ref_ptr<Channels> channels) {
    this->channels = channels;
}

void NeutrinoAntineutrinoInteraction::setChannelsBundle(ref_ptr<Channels> channels, std::string fname) {
    this->channelsBundle = new ChannelsBundle(channels, fname);
}

void NeutrinoAntineutrinoInteraction::setRelativisticInteraction(double m1, double m2, double E, double s) const {
    this->relInteraction = new RelativisticInteraction(m1, m2, E, s);
}

// the structure is inspired by EMInverseCompton class etc.
class NeutrinoAntineutrinoSecondariesDistribution {
    private:
        std::vector< std::vector<double> > data;
        std::vector<double> s_values;
        std::vector<double> costh13_com_values;
        size_t Ns;
        size_t Nrer;
        double s_min;
        double s_max;
        double costh_min;
        double costh_max;
        double dcosth;
        double dls;
    
    public:
        
    NeutrinoAntineutrinoSecondariesDistribution(std::string variable, int idChannel, double sThr) {
        
        Ns = 1000;
        Nrer = 1000;
        
        s_min = sThr; // yet in eV
        s_max = 1e28 * eV * eV; // value from the computed tables
        dls = (log(s_max) - log(s_min)) / Ns;
        
        costh_min = -1;
        costh_max = 1;
        dcosth = (costh_max - costh_min) / Nrer;
        
        data = std::vector< std::vector<double> >(1000, std::vector<double>(1000));
        std::vector<double> data_i(1000);
        
        // tabulating the s bin borders
        s_values = std::vector<double>(1001);
        for (size_t i = 0; i < Ns + 1; ++i)
            s_values[i] = s_min * exp(i * dls);
        
        if (variable != "costh13_com")
            throw std::runtime_error("The only available variable to compute the differential cross section is costheta13_com!");
        
        // tabulating the costh13_com bin borders
        costh13_com_values = std::vector<double>(1001);
        for (size_t i = 0; i < Nrer + 1; ++i)
            costh13_com_values[i] = costh_min + i * dcosth;
        
        // to generate randomly the seed to produce the differential cross section
        // Random &random = Random::instance();
        // int seedDiffXS = random.randInt();
        // not affecting the results
        int seedDiffXS = 1;
        
        // for each s and costh13_com tabulate the cumulative differential cross section
        for (size_t i = 0; i < Ns; i++) {
            double s = s_min * exp((i + 0.5) * dls);
            
            // cumulative midpoint integration
            data_i[0] = getDifferentialXS(s, variable, costh_min, idChannel, seedDiffXS) * dcosth * i;
            // to check * (dcosth[i] - 1)!!!
            
            for (size_t j = 1; j < Nrer; j++) {
                double costh13_com = costh_min + (j + 0.5) * dcosth;
                double dcosth13_com = (j + 1) * dcosth - j * dcosth; // since we are on a linear scale, dcosth13_com == dcosth
                
                data_i[j] = getDifferentialXS(s, variable, costh13_com, idChannel, seedDiffXS) * dcosth13_com;
                data_i[j] += data_i[j - 1];
            }
            data[i] = data_i;
        }
        
    }
    
    double getDifferentialXS(double s, std::string variable, double variableValue, int idChannel, int seedDiffXS) {
        
        std::string partonicPath = "/Applications/CRPropa/NuPropaLap/PartonicCalculation/sigmaNu_interface/"; // to change with NUPROPA path
        std::string interfacePath = partonicPath + "bin/";

        std::ostringstream cmd;
        cmd << interfacePath << "Main_Interface.exe"
            << " -c " << idChannel
            << " -s " << seedDiffXS
            << " -E " << std::sqrt(s / GeV / GeV); // Ecms has to be given in GeV
        
        int result = std::system(cmd.str().c_str());
        if (result != 0) {
            std::ostringstream oss;
            oss << "Error: Failed to run Main_Interface.exe, with exit code: " << result;
            throw std::runtime_error(oss.str());
        }

        std::ostringstream ss;
        ss << std::scientific << std::setprecision(5) << s / GeV / GeV;
        std::string Ecms2 = ss.str();
        
        std::string filePath = partonicPath + "dataDifferentialXS/channel" + std::to_string(idChannel) + "/";
        std::string filename = filePath + variable + "_channel" + std::to_string(idChannel) + "_EcmsSq" + Ecms2 + "_s" + std::to_string(seedDiffXS) + ".txt";
        
        std::vector<double> variableScan;
        std::vector<double> differentialXS;
        
        std::ifstream infile(filename);
        
        if(!infile) {
            std::ostringstream oss;
            oss << "Error: could not open the differentialXS file. The filename is: " << filename;
            throw std::runtime_error(oss.str());
        } else {
            
            double a, b;
            while (infile >> a >> b) {
                variableScan.push_back(a);
                differentialXS.push_back(b);
            }
            infile.close();
        }
        
        double diffXS = interpolate(variableValue, variableScan, differentialXS);
        return diffXS;
        
    }
    
    // draw random costh13_com from the differential cross section distribution
    double sample(double s) {
        size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
        std::vector<double> s0 = data[idx];
        
        Random &random = Random::instance();
        size_t j = random.randBin(s0); // draw random bin (lower bin boundary returned)
        double costh13_comRand = costh_min + j * dcosth;
        
        return costh13_comRand;
    }
        
};

void NeutrinoAntineutrinoInteraction::performInteraction(Candidate *candidate, double mass, int IDbkg) const {
    
    double E = candidate->current.getEnergy();
    int ID = candidate->current.getId();
    
    candidate->setActive(false);
    
    if (not haveSecondaries)
        return;
    
    int indexChannel = this->channelsBundle->getSelectedIndex();
    
    std::vector<double> tabE = this->channelsBundle->selectE();
    std::vector<double> tabs = this->channelsBundle->selects();
    std::vector<std::vector<double>> tabCDF = this->channelsBundle->selectCDF();
    std::vector<int> prodChanId = this->channelsBundle->selectProdChanId(); // the channel ID is inside the products.txt file in each folder.
    
    // check if in tabulated energy range
    if (E < tabE.front() or (E > tabE.back()))
        return;
    
    // solve the problem for elastic scattering, in this case I should take ID & IDbkg
    int prodId1 = prodChanId[0];
    int prodId2 = prodChanId[1];
    int chanId = prodChanId[2];
    
    double m3;
    double m4;
    
    if (prodId1 == 0 && prodId2 == 0) {  // in the elastic scattering they are set to 0!
        
        prodId1 = prodId1 + ID;
        prodId2 = prodId2 + IDbkg;
        
        m3 = mass;
        m4 = this->neutrinoFieldMass;
        
    } else {
        
        ParticleData particle;
        m3 = particle.getParticleMass(prodId1); // in kg
        m4 = particle.getParticleMass(prodId2); // in kg
        
    }
    
    double sThr = (m3 * m3 + m4 * m4) * c_squared * c_squared;
    
    // sample the value of s
    Random &random = Random::instance();
    size_t i = closestIndex(E, tabE);  // find closest tabulation point
    size_t j = random.randBin(tabCDF[i]);
    double lo = std::max(sThr, tabs[j-1]); // first s-tabulation point below min(s_kin);
    double hi = tabs[j];
    double s = lo + random.rand() * (hi - lo); // should I add the neutrino masses? since it is the s_kin!!
    
    // sample the cosine of theta13_com
    std::string variable = "costh13_com";
    static NeutrinoAntineutrinoSecondariesDistribution distribution(variable, chanId, sThr);
    double costh13_com = distribution.sample(s);

    setRelativisticInteraction(mass, this->neutrinoFieldMass, E, s); // to treat for the case one or both massless
    
    // energies of the secondary particles
    std::vector<double> energies = this->relInteraction->getProductEnergiesLab(s, costh13_com, m3, m4);
    
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    
    // add the thinning!
    double w = 1.;
    
    candidate->addSecondary(prodId1, energies[0], pos, w, this->interactionTag);
    candidate->addSecondary(prodId2, energies[1], pos, w, this->interactionTag);
    
}

void NeutrinoAntineutrinoInteraction::process(Candidate *candidate) const {
    
    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy();
    int ID = candidate->current.getId();
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;

    int IDbkg = this->neutrinoMixing->fromMassToFlavour(this->neutrinoFieldMass / eV);
    
    if (ID * IDbkg > 0)
        return;
    
    double mass = this->neutrinoMixing->fromFlavourToMass(ID) * eV; // returned in eV from the function
    int indexMass = this->neutrinoMixing->massToIndexMass(mass / eV) + 1;
    std::string massComb = this->neutrinoField->getFieldName() + "_m" + std::to_string(indexMass);
    
    double rate = this->channelsBundle->getRate(ID, IDbkg, massComb, z, E);
    
    if (rate < 0)
        return;
    
    // check for interaction
    Random &random = Random::instance();
    double randDistance = -log(random.rand()) / rate;
    double step = candidate->getCurrentStep();
    if (step < randDistance) {
        candidate->limitNextStep(limit / rate);
        return;
    } else { // after performing interaction neutrino ceases to exist (hence return)
        performInteraction(candidate, mass, IDbkg);
        return;
    }
}

void NeutrinoAntineutrinoInteraction::setInteractionTag(std::string tag) const {
    this->interactionTag = tag;
}

std::string NeutrinoAntineutrinoInteraction::getInteractionTag() const {
    return this->interactionTag;
}

}

