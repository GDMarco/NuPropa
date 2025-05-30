#include "nupropa/NeutrinoPhotonInteraction.h"
#include "nupropa/NeutrinoBackground.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Module.h>
#include <crpropa/Candidate.h>

#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <filesystem>

namespace nupropa {

using namespace crpropa;

// The parent's constructor need to be called on initialization!
NeutrinoPhotonInteraction::NeutrinoPhotonInteraction(ref_ptr<PhotonField> photonField, bool haveSecondaries,  double limit) : Module() { //double thinning,
    setPhotonField(photonField);
    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    //setThinning(thinning);
    
}

void NeutrinoPhotonInteraction::setPhotonField(ref_ptr<PhotonField> photonField) {
    this->photonField = photonField;
    std::string fname = photonField->getFieldName();
    setDescription("NeutrinoPhotonInteraction::Module" + fname);
    initRate(fname);
}

void NeutrinoPhotonInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoPhotonInteraction::setLimit(double limit) {
    this->limit = limit;
}

/**
void NeutrinoPhotonInteraction::setThinning(double thinning) {
    this->thinning = thinning;
}
*/

void NeutrinoPhotonInteraction::initRate(std::string fname) {
    
    tabEnergy.clear();
    tabRate.clear();
    
    std::unordered_map<int, std::string> dictionaryNeutrino;
    std::__fs::filesystem::path dir = "/Applications/CRPropa/NuGammaInteraction/CRPropa3-data/data/NeutrinoInteractions/NeutrinoPhotonInteraction/";
    //getDataPath(getDataPath("") + "data/NeutrinoPhotonInteraction/"); // to check!
    
    int index = 0;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {
        
        std::string filePath = dir_entry.path().string() + "/rate_" + fname + ".txt"; // to check!
        std::ifstream infile(filePath.c_str());
        
        if (!infile.good())
            throw std::runtime_error("NeutrinoPhotonInteraction: could not open file" + filePath);
        
        std::vector<double> vecEnergy;
        std::vector<double> vecRate;
        
        while (infile.good()) {
            if (infile.peek() != '#') {
                double a, b;
                infile >> a >> b;
                if (infile) {
                    vecEnergy.push_back(pow(10, a) * eV);
                    vecRate.push_back(b / Mpc);
                }
            }
            infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
        }
        infile.close();
        
        this->tabEnergy.push_back(vecEnergy);
        this->tabRate.push_back(vecRate);

        dictionaryNeutrino[index] = dir_entry.path().filename().string();
        index = index + 1;
    }
    this->dictionaryNeutrino = dictionaryNeutrino;
}

void NeutrinoPhotonInteraction::performInteraction(Candidate *candidate) const {
    
    candidate->setActive(false);

    if (not haveSecondaries)
        return;

    double w = 1.; // no weights for now!
    // Use assumption of Seckel, 1998: W boson production on-shell.
    // W boson produced on shell, mass_W needs to be defined in Units.h:
    // static const double mass_W = 80.377 * 1e9 * 1.602176487e-19 *  (2.99792458e-2 * 1e-16) * kilogram;
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    int ID = candidate->current.getId();
    
    double leptonE;
    
    if (std::abs(ID) == 12) {
        leptonE = (E - (mass_W + mass_electron) * c_squared); // to adjust according to the ID!
    } else if (std::abs(ID) == 14) {
        leptonE = (E - (mass_W + mass_muon) * c_squared);
    } else {
        leptonE = (E - (mass_W + mass_tauon) * c_squared);
    }
    
    int leptonID;
    
    if (ID > 0) {
        leptonID = ID - 1;
    } else {
        leptonID = ID + 1;
    }

    Random &random = Random::instance();
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    
    if (haveSecondaries)
        candidate->addSecondary(leptonID, leptonE / (1 + z), pos, w, interactionTag);
        //candidate->addSecondary(W)
}

std::vector<double> NeutrinoPhotonInteraction::getTabulatedEnergy(int ID) const {
    
    int indexInteraction;

    if (abs(ID) == 12) {
        std::string interaction = "NeutrinoElectronPhotonInteraction";
        for (const auto& pair : this->dictionaryNeutrino) {
            if (pair.second == interaction) {
                    indexInteraction = pair.first;
                    break;
                }
        }
    } else if (abs(ID) == 14) {
        std::string interaction = "NeutrinoMuonPhotonInteraction";
        for (const auto& pair : this->dictionaryNeutrino) {
                if (pair.second == interaction) {
                    indexInteraction = pair.first;
                    break;
                }
        }
    } else if (abs(ID) == 16) {
        std::string interaction = "NeutrinoTauPhotonInteraction";
        for (const auto& pair : this->dictionaryNeutrino) {
            if (pair.second == interaction) {
                indexInteraction = pair.first;
                break;
            }
        }
    } else {
        throw std::runtime_error("Invalid ID: tables not found!");
    }
    return this->tabEnergy[indexInteraction];
}

std::vector<double> NeutrinoPhotonInteraction::getTabulatedRate(int ID) const {
    
    int indexInteraction; // think about initialising to a random value, i.e. 5
    
    if (abs(ID) == 12) {
        std::string interaction = "NeutrinoElectronPhotonInteraction";
        for (const auto& pair : this->dictionaryNeutrino) {
                if (pair.second == interaction) {
                    indexInteraction = pair.first;
                    break;
                }
            }
    } else if (abs(ID) == 14) {
        std::string interaction = "NeutrinoMuonPhotonInteraction";
        for (const auto& pair : this->dictionaryNeutrino) {
                if (pair.second == interaction) {
                    indexInteraction = pair.first;
                    break;
                }
            }
    } else if (abs(ID) == 16) {
        std::string interaction = "NeutrinoTauPhotonInteraction";
        for (const auto& pair : this->dictionaryNeutrino) {
            if (pair.second == interaction) {
                indexInteraction = pair.first;
                break;
            }
        }
    } else {
        throw std::runtime_error("Invalid ID: tables not found!");
    }
    return this->tabRate[indexInteraction];
}

void NeutrinoPhotonInteraction::process(Candidate *candidate) const
{
    // To enable parallelization, the modules have to be stateless - the
    // process method should thus not modify internal variables!
    //std::cout << "NeutrinoPhotonInteraction::Module::process() called\n";
    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = (1 + z) * candidate->current.getEnergy();
    double ID = candidate->current.getId();
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;
   
    std::vector<double> tabEnergy = getTabulatedEnergy(ID);
    std::vector<double> tabRate = getTabulatedRate(ID);
    
    // check if in tabulated energy range
    if (E < tabEnergy.front() or (E > tabEnergy.back()))
        return;

    // interaction rate
    double rate = interpolate(E, tabEnergy, tabRate);
    rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

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

void NeutrinoPhotonInteraction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string NeutrinoPhotonInteraction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace nupropa
