#include "NeutrinoAntineutrinoInteraction.h"
#include "Channels.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <string>
#include <filesystem>
#include <unordered_map>

namespace crpropa {

// The parent's constructor need to be called on initialization!
NeutrinoAntineutrinoInteraction::NeutrinoAntineutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, bool haveSecondaries, ref_ptr<Channels> channels, double limit) : Module() { // double thinning
    setChannels(channels);
    setNeutrinoField(neutrinoField);
    setHaveSecondary(haveSecondaries);
    setLimit(limit);
    //setThinning(thinning);
    // 6 channels for the two interactions
}

/**
 One assumes that the tables have been produced within the same computation, so the energy bins are the same.
 
 In the case of the neutrino-antineutrino interactions the elastic scatterings furnish the tabEnergy, since they do not have interaction energy thresholds!
 */

// channels should be an unordered map <bool, "interaction">, maybe interaction tag? setChannels should be a way of filling this object, knowing how they are ordered!
void NeutrinoAntineutrinoInteraction::setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField) {
    this->neutrinoField = neutrinoField;
    this->neutrinoFieldID = neutrinoField->getParticleID(); 
    std::string fname = neutrinoField->getFieldName();
    
    setDescription("NeutrinoAntineutrinoInteraction::Module" + fname);
    
    initRate(fname);
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

void NeutrinoAntineutrinoInteraction::initRate(std::string fname) {
    
    tabEnergy.clear();
    tabRate.clear();
    
    std::string pathInteraction = this->channels->getInteractionFolderPath();
    std::vector<std::string> interactionChannel = this->channels->getInteractionChannels();
    std::vector<bool> activeChannels = this->channels->getActiveChannels();
    std::unordered_map<int, std::string> interactionDictionary;
    int iChannels = 0;
    
    for (int i; i <= interactionChannel.size(); i++) {
            
        if (!activeChannels[i]) {
            continue;
        } else {
            filename = path + interactionChannels[i] + "/rate_" + fname + ".txt";
            std::ifstream infile(filename.c_str());
            
            
            
            if (!infile.good())
                throw std::runtime_error(folder + ": could not open file" + filename);
            
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
            interactionDictionary[iChannels] = interactionChannels[i];
            
            iChannels = iChannels + 1;
        }
    
        this->interactionDictionary = interactionDictionary;
    }
}

int NeutrinoAntineutrinoInteraction::searchChannel(std::string interacting, std::string products) {
    std::string totalChannel = interacting + products;
    
    for (const auto& el : this->interactionDictionary) {
        std::string interactionChannel = el.second;
        
        if (interactionChannel == totalChannel) {
            return el.first;
            break;
        } else {
            continue;
        }
    }
    return 50; // not active channel
}

std::vector<double> NeutrinoAntineutrinoInteraction::fillTableZeros(std::vector<double> table, size_t size) {
    if (table.size() < size) {
        difference = size - table.size();
        table.insert(table.begin(), difference, 0);
        return resizedTable;
    } else if (table.size() == size) {
        return table;
    } else {
        throw std::runtime_error("Too large table size: suspicious!");
    }
}

std::pair<std::vector<double>, std::vector<double>> NeutrinoAntineutrinoInteraction::getRateTables(int ID, int nuBkgID) {
    if (abs(ID) == abs(nuBkgID)) {
        setInteractionTag("NuAntiNuInt");
        
        int iElastic = searchChannel("NeutrinoAntineutrino", "Elastic"); // negative values in the tables! TO CHECK!
        int iWProd = searchChannel("NeutrinoAntineutrino", "WProduction");
        int iLep;
        int iResEl = searchChannel("NeutrinoAntineutrino", "ResonanceElectron");
        int iResMu = searchChannel("NeutrinoAntineutrino", "ResonanceMuon");
        int iResTa = searchChannel("NeutrinoAntineutrino", "ResonanceTau");
        int iResUp = searchChannel("NeutrinoAntineutrino", "ResonanceUp");
        int iResDown = searchChannel("NeutrinoAntineutrino", "ResonanceDown");
        int iResCharm = searchChannel("NeutrinoAntineutrino", "ResonanceCharm");
        int iResStrange = searchChannel("NeutrinoAntineutrino", "ResonanceStrange");
        int iResTop = searchChannel("NeutrinoAntineutrino", "ResonanceTop");
        int iResBottom = searchChannel("NeutrinoAntineutrino", "ResonanceBottom");
        
        if (abs(ID) == 12) {
            iLep = searchChannel("NeutrinoAntineutrino", "Electron");
        } else if (abs(ID) == 14) {
            iLep = searchChannel("NeutrinoAntineutrino", "Muon");
        } else {
            iLep = searchChannel("NeutrinoAntineutrino", "Tau");
        }
        
        std::vector<int> indexes = {iElastic, iWProd, iLep, iResEl, iResMu, iResTa, iResUp, iResDown, iResCharm, iResStrange, iResTop, iResBottom}; // fixed vector for all the possible channels in NeutrinoAntineutrinoInteraction
        size_t maximumTableSize = 0;
        int maximumIndex;
        std::vector<int> activeIndexes;
        
        for (size_t j = 0; j < indexes.size(); j++) {
            if (indexes[j] == 50) {
                continue;
            } else {
                int tableSize = this->tabRate[indexes[j]];
                if (tableSize > maximumTableSize) {
                    maximumTableSize = tableSize;
                    maximumIndex = indexes[j];
                }
                activeIndexes.push_back(indexes[j]);
            }
        }
        
        std::vector<double> totalRate(maximumTableSize, 0);
        
        for (size_t j; j < activeIndexes.size(); j++) {
            std::vector<double> resizedRate = fillTableZeros(this->tabRate[activeIndexes[j]])
            
            for (int k = 0; k < maximumTableSize; k++) {
                totalRate[k] = totalRate[k] + resizedRate[k];
            }
        }
        
        return std::make_pair<this->tabEnergy[maximumIndex], totalRate>;
        
    } else {
        setInteractionTag("NuiAntiNujInt");
        
        int iElastic = searchChannel("NeutrinoiAntineutrinoj", "Elastic");
        std::vector<double> tabRateElastic = this->tabRate[iElastic];
        
        int iLep;
        
        int totalFlavour = abs(ID) + abs(nuBkgID);
        if (totalFlavour == 26) { // for the rate there's no difference on the cross section apart from the energy threshold that depends on the mass of the products, TO CHECK!
            iLep = searchChannel("NeutrinoiAntineutrinoj", "ElectronAntimuon");
        } else if (totalFlavour == 28) {
            iLep = searchChannel("NeutrinoiAntineutrinoj", "ElectronAntitau");
        } else {    // i.e. else if (totalFlavour == 30)
            iLep = searchChannel("NeutrinoiAntineutrinoj", "TauAntimuon");
        }
        
        if (iElastic == 50) {
            return std::make_pair<this->tabEnergy[iLep], this->tabRate[iLep]>;
        } else {
            size_t sizeElasticRate = this->tabRate[iElastic].size();
            std::vector<double> elasticRate = this->tabRate[iElastic];
            std::vector<double> resizedLepRate = fillTableZeros(this->tabRate[iLep], sizeElasticRate);
            
            std::vector<double> totalRate(sizeElasticRate, 0);
            for (int k = 0; k < sizeElasticRate; k++) {
                totalRate[k] = elasticRate[k] + resizedLepRate[k];
            }
            return std::make_pair<this->tableEnergy[iElastic], totalRate>
        }
    }
}

void NeutrinoAntineutrinoInteraction::performInteraction(Candidate *candidate) const {
    
    candidate->setActive(false);
    
    if (not haveSecondaries)
        return;
    
    // function to take the cumulative inverse mean free path and the other stuff!
    
    // the products depends on the channel
    double w = 1.;
    
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    int ID = candidate->current.getID();
    
    // it is better to take the ID of the particles produced from the Channels class!
    // IDs of the secondaries particles
    // energies of the secondary particles
    Random &random = Random::instance();
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    
    if (haveSecondaries)
        candidate->addSecondary(leptonID, leptonE / (1 + z), pos, w, interactionTag);
        //candidate->addSecondary()
}

void NeutrinoAntineutrinoInteraction::process(Candidate *candidate) const
{
    
    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = (1 + z) * candidate->current.getEnergy();
    int ID = candidate->current.getID();
    int nuBkgID = this->neutrinoFieldID;
    
    if (!(abs(ID) == 12 || abs(ID) == 12 || abs(ID) == 16))
        return;

    auto tablesRate = getRateTables(ID, nuBkgID); // another option is to use a struct twoTables...
    std::vector<double> tabEnergy = tablesRate.first;
    std::vector<double> tabRate = tablesRate.second;
    
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

void NeutrinoAntineutrinoInteraction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string NeutrinoAntineutrinoInteraction::getInteractionTag() const {
    return interactionTag;
}

}

