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
    
    for (int i; i <= interactionChannel.size(); ++i) {
            
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
            this->tabProductsID.push_back(this->channels->getProductsID()[i]) // see if it works, if not define it in advance
            
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

// think about moving it to Channels.cc, using as arguments this->tabEnergy(), this->tabRate()
std::pair<std::vector<double>, std::vector<double>> NeutrinoAntineutrinoInteraction::getRateTables(int ID, int IDBkg) {
    if (abs(ID) == abs(IDBkg)) {
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
                this->tabProductsID.erase(this->tabProductsID.begin() + j);
            } else {
                // maybe also here is needed the erase function! more efficient
                int tableSize = this->tabRate[indexes[j]];
                if (tableSize > maximumTableSize) {
                    maximumTableSize = tableSize;
                    maximumIndex = indexes[j];
                }
                activeIndexes.push_back(indexes[j]);
            }
        }
        
        std::vector<double> totalRate(maximumTableSize, 0);
        std::vector<std::vector<double>> resizedVectors;
        
        for (size_t j; j < activeIndexes.size(); ++j) {
            std::vector<double> resizedRate = fillTableZeros(this->tabRate[activeIndexes[j]], maximumTableSize)
            resizedVectors.append(resizedRate);
            
            for (int k = 0; k < maximumTableSize; ++k) {
                totalRate[k] = totalRate[k] + resizedRate[k];
            }
        } // to review!
        
        computeInteractionProbabilities(resizedVectors); // make sure it is sync with the interactionDictionary --> it's not sync with the interaction dictionary.
        return std::make_pair<this->tabEnergy[maximumIndex], totalRate>;
        
    } else {
        
        setInteractionTag("NuiAntiNujInt");
        
        int iElastic = searchChannel("NeutrinoiAntineutrinoj", "Elastic");
        std::vector<double> tabRateElastic = this->tabRate[iElastic];
        
        int iLep;
        
        int totalFlavour = abs(ID) + abs(IDBkg);
        if (totalFlavour == 26) { // for the rate there's no difference on the cross section apart from the energy threshold that depends on the mass of the products, TO CHECK!
            // now it does matter! to establish for the six cases! parity equivalence
            // actually it does not matter, since to find the products I use another method
            iLep = searchChannel("NeutrinoiAntineutrinoj", "ElectronAntimuon");
        } else if (totalFlavour == 28) {
            iLep = searchChannel("NeutrinoiAntineutrinoj", "ElectronAntitau");
        } else {    // i.e. else if (totalFlavour == 30)
            iLep = searchChannel("NeutrinoiAntineutrinoj", "TauAntimuon");
        }
        
        if ((iElastic == 50) && (iLep != 50)) { // the case for iLep == 50??? and for both 50??
            size_t size = this->tabEnergy[iLep].size();
            std::vector<double> vecProb(size, 1);
            this->channelProbability.push_back(vecProb);
            
            return std::make_pair<this->tabEnergy[iLep], this->tabRate[iLep]>;
        } else if ((iElastic != 50) && (iLep == 50)) {
            size_t size = this->tabEnergy[iElastic].size();
            std::vector<double> vecProb(size, 1);
            this->channelProbability.push_back(vecProb);
            
            return std::make_pair<this->tabEnergy[iElastic], this->tabRate[iElastic]>;
            
        } else if ((iElastic != 50) && (iLep != 50)) {
            
            size_t sizeElasticRate = this->tabRate[iElastic].size();
            std::vector<double> elasticRate = this->tabRate[iElastic];
            std::vector<double> resizedLepRate = fillTableZeros(this->tabRate[iLep], sizeElasticRate);
            
            std::vector<double> totalRate(sizeElasticRate, 0);
            std::vector<std::vector<double>> resizedVectors;
            
            for (int k = 0; k < sizeElasticRate; k++) {
                totalRate[k] = elasticRate[k] + resizedLepRate[k];
            }
            
            resizedVectors.push_back(resizedLepRate);
            resizedVectors.push_back(elasticRate);
            computeInteractionProbability(resizedVectors);
            
            return std::make_pair<this->tableEnergy[iElastic], totalRate>;
            
        } else {
            throw std::runtime_error("No active nu-antinu channels! Reset it.");
        }
    }
}

void NeutrinoAntineutrinoInteraction::computeInteractionProbabilities (std::vector<std::vector<double>> rates) {
    size_t cols = rates[0].size();
    size_t rows = rates.size();
    
    for (size_t i; i <= cols; ++i) {
        double sum = 0;
        for (size_t j; j <= rows; ++j) {
            sum = sum + this->tabRate[j][i];
        }
        
        std::vector<double> vecProb;
        
        for (size_t j; j <= rows; ++j) {
            vecProb.push_back(this->tabRate[j][i] / sum);
        }
        
        this->channelProbability.push_back(vecProb);
        // to check if it's work, and to sync with interaction dictionary!
    }
    // in the file it should read I need both the products and the interactive particles, to
    //return this->channelProbabilities();
    // tabProbabilities // channels (index) vs pbin1 - pbin1 - pbin2
                        // 1                    0.1
                        // 2                    0.8
                        // 3                    0.1
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
    
    // IDs of the secondaries particles
    std::vector<int> getProductsID(ID, IDBkg);
    
    // energies of the secondary particles
    Random &random = Random::instance();
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    
    if (haveSecondaries)
        candidate->addSecondary(ID1, E1 / (1 + z), pos, w, interactionTag);
        candidate->addSecondary(ID2, E2 / (1 + z), pos, w, interactionTag);
}

void NeutrinoAntineutrinoInteraction::process(Candidate *candidate) const
{
    
    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = (1 + z) * candidate->current.getEnergy();
    int ID = candidate->current.getID();
    int IDBkg = this->neutrinoFieldID;
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;

    auto tablesRate = getRateTables(ID, IDBkg); // another option is to use a struct twoTables...
    std::vector<double> tabEnergy = tablesRate.first;
    std::vector<double> tabRate = tablesRate.second;
   
    // check if in tabulated energy range
    if (E < tabEnergy.front() or (E > tabEnergy.back()))
        return;
    
    //from here it might be a function to get the channel!
    // take the index i of the closest tabEnergy[i] to E
    double distEnergy = 1e6; // to minimize this value
    int iEnergy = -1;
    for (int i; i <= tabEnergy.size(); ++i) {
        dist = abs(E - tabEnergy[i]);
        if (dist < distEnergy) {
            distEnergy = dist;
            iEnergy = i;
        }
    }
    
    // select the channel and give the IDs
    std::vector<double> channelProbability = this->channelProbability[iEnergy];
    
    // randomly select the interaction channel
    Random &randomChannel = Random::instance();
    double channel = randomChannel.rand();
    
    double tabCDFLow = 0;
    double tabCDFUpp = 0;
    int selChannel = -1;
    
    // build a "cumulative function" and select the channel
    for (int i; i <= channelProbability.size() - 1; ++i) {
        tabCDFLow = tabCDFLow + channelProbability[i];
        tabCDFUpp = tabCDFUpp + channelProbability[i + 1];
        if (tabCDFLow < iRand <= tabCDFUpp) {
            selChannel = i;
            break;
        }
    }
    
    std::vector<int> IDs = this->tabProductsID[selChannel];
    // see the correspondence between the id and the selected channel -> assign the products IDs
    // to take into account parity invariance of the processes... !
    // for the elastic scatterings there are initialising values!
    
    // end of the function to get the IDs!
    
    // interaction rate
    double rate = interpolate(E, tabEnergy, tabRate);
    rate *= pow_integer<2>(1 + z) * neutrinoField->getRedshiftScaling(z); // see how it is implemented this function in NeutrinoBackground.h
    
    
    
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

