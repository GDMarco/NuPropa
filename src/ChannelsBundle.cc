#include "nupropa/ChannelsBundle.h"
#include "nupropa/Channels.h"

#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Random.h>
#include <crpropa/Common.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>

namespace nupropa {

using namespace crpropa;

ChannelsBundle::ChannelsBundle() {};

ChannelsBundle::ChannelsBundle(ref_ptr<Channels> channels, std::string fname) {
    
    std::string path = channels->getInteractionFolderPath();
    std::vector<std::vector<int>> products = channels->getProductsID();
    std::vector<std::string> interactionChannel = channels->getInteractionChannels();
    std::vector<bool> activeChannels = channels->getActiveChannels();
    
    std::unordered_map<int, std::string> interactionDictionary;
    int iChannels = 0;
    
    for (int i = 0; i <= interactionChannel.size(); ++i) {
        if (!activeChannels[i]) {
            continue;
        } else {
            std::string filename = path + interactionChannel[i] + "/rate_" + fname + ".txt";
            std::ifstream infile(filename.c_str());
            
            if (!infile.good())
                throw std::runtime_error("Could not open rate file: " + filename);
            
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
            
            std::string filenameCDF = path + interactionChannel[i] + "/cdf_" + fname + ".txt";
            
            std::ifstream infileCDF(filenameCDF.c_str());
            
            if (!infileCDF.good())
                throw std::runtime_error("Could not open CDF file: " + filenameCDF);
            
            std::vector<double> vecE;
            std::vector<double> vecs;
            std::vector<std::vector<double>> vecCDF;
            
            // skip header
            while (infileCDF.peek() == '#')
                infileCDF.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
            
            double c;
            // read s values in first line
            infileCDF >> c; // skip first value
            while (infileCDF.good() and (infileCDF.peek() != '\n')) {
                infileCDF >> c;
                vecs.push_back(pow(10, c) * eV * eV);
            }
            
            // read all following lines: E, cdf values
            while (infileCDF.good()) {
                infileCDF >> c;
                if (!infileCDF)
                    break;  // end of file
                vecE.push_back(pow(10, c) * eV);
                std::vector<double> cdf;
                for (int i = 0; i < vecs.size(); i++) {
                    infileCDF >> c;
                    cdf.push_back(c / Mpc);
                }
                vecCDF.push_back(cdf);
            }
            infileCDF.close();
            
            this->tabEnergy.push_back(vecEnergy);
            this->tabRate.push_back(vecRate);
            
            this->tabE.push_back(vecE);
            this->tabs.push_back(vecs);
            this->tabCDF.push_back(vecCDF);
            
            this->tabProductsID.push_back(products[i]); // see if it works, if not define it in advance
            
            interactionDictionary[iChannels] = interactionChannel[i];
            
            iChannels = iChannels + 1;
        }
        this->interactionDictionary = interactionDictionary;
    }
}

void ChannelsBundle::sortDictionaryIndexes(int indexErased) {
    
    std::unordered_map<int, std::string> tempDictionary;
    
    for (auto el = this->interactionDictionary.begin(); el != this->interactionDictionary.end(); ++el) {
        
        int currentIndex = el->first;
        std::string currentInteraction = el->second;
        
        if (currentIndex == indexErased) {
            continue;
        } else if (el->first > indexErased) {
            tempDictionary[currentIndex - 1] = currentInteraction;
        } else {
            tempDictionary[currentIndex] = currentInteraction;
        }
    }
    
    this->interactionDictionary = tempDictionary;
    
}

void ChannelsBundle::removeChannels(int ID, int IDBkg) {
    
    if (abs(ID) == abs(IDBkg)) {
        
        std::string interacting = "NeutrinoAntineutrino";
        
        for (const auto& el : this->interactionDictionary) {
            std::string interactionChannel = el.second;
            
            if (interactionChannel.find(interacting) == 0) {
                
                int i = el.first;
                
                if (abs(ID) == 12) {
                    std::string leptonic1 = interacting + "Muon";
                    std::string leptonic2 = interacting + "Tau";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                    
                        sortDictionaryIndexes(i);
                    }
                    
                } else if (abs(ID) == 14) {
                    std::string leptonic1 = interacting + "Electron";
                    std::string leptonic2 = interacting + "Tau";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                } else {
                    std::string leptonic1 = interacting + "Muon";
                    std::string leptonic2 = interacting + "Electron";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                }
            } else {
                
                int i = el.first;
                
                this->tabEnergy.erase(this->tabEnergy.begin() + i);
                this->tabRate.erase(this->tabRate.begin() + i);
                this->tabE.erase(this->tabE.begin() + i);
                this->tabs.erase(this->tabs.begin() + i);
                this->tabCDF.erase(this->tabCDF.begin() + i);
                this->tabProductsID.erase(this->tabProductsID.begin() + i);
                
                sortDictionaryIndexes(i);
            }
        }
   
    } else {
        
        std::string interacting = "NeutrinoiAntineutrinoj";
    
        for (const auto& el : this->interactionDictionary) {
            std::string interactionChannel = el.second;
            
            if (interactionChannel.find(interacting) == 0) {
                
                int i = el.first;
                
                if (abs(ID) == 12 && abs(IDBkg) == 14) {
                    std::string leptonic1 = interacting + "MuonAntielectron";
                    std::string leptonic2 = interacting + "TauAntielectron";
                    std::string leptonic3 = interacting + "ElectronAntitau";
                    std::string leptonic4 = interacting + "MuonAntitau";
                    std::string leptonic5 = interacting + "TauAntimuon";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2 || interactionChannel == leptonic3 || interactionChannel == leptonic4 || interactionChannel == leptonic5) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                } else if (abs(ID) == 12 && abs(IDBkg) == 16) {
                    std::string leptonic1 = interacting + "MuonAntielectron";
                    std::string leptonic2 = interacting + "TauAntielectron";
                    std::string leptonic3 = interacting + "ElectronAntimuon";
                    std::string leptonic4 = interacting + "MuonAntitau";
                    std::string leptonic5 = interacting + "TauAntimuon";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2 || interactionChannel == leptonic3 || interactionChannel == leptonic4 || interactionChannel == leptonic5) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                } else if (abs(ID) == 14 && abs(IDBkg) == 12) {
                    std::string leptonic1 = interacting + "ElectronAntimuon";
                    std::string leptonic2 = interacting + "TauAntielectron";
                    std::string leptonic3 = interacting + "ElectronAntimuon";
                    std::string leptonic4 = interacting + "MuonAntitau";
                    std::string leptonic5 = interacting + "TauAntimuon";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2 || interactionChannel == leptonic3 || interactionChannel == leptonic4 || interactionChannel == leptonic5) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                } else if (abs(ID) == 14 && abs(IDBkg) == 16) {
                    std::string leptonic1 = interacting + "ElectronAntimuon";
                    std::string leptonic2 = interacting + "TauAntielectron";
                    std::string leptonic3 = interacting + "ElectronAntimuon";
                    std::string leptonic4 = interacting + "MuonAntielectron";
                    std::string leptonic5 = interacting + "TauAntimuon";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2 || interactionChannel == leptonic3 || interactionChannel == leptonic4 || interactionChannel == leptonic5) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                } else if (abs(ID) == 16 && abs(IDBkg) == 12) {
                    std::string leptonic1 = interacting + "ElectronAntimuon";
                    std::string leptonic2 = interacting + "MuonAntitau";
                    std::string leptonic3 = interacting + "ElectronAntimuon";
                    std::string leptonic4 = interacting + "MuonAntielectron";
                    std::string leptonic5 = interacting + "TauAntimuon";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2 || interactionChannel == leptonic3 || interactionChannel == leptonic4 || interactionChannel == leptonic5) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                } else {
                    // abs(ID) == 16 && abs(IDBkg) == 14
                    std::string leptonic1 = interacting + "ElectronAntimuon";
                    std::string leptonic2 = interacting + "TauAntielectron";
                    std::string leptonic3 = interacting + "ElectronAntimuon";
                    std::string leptonic4 = interacting + "MuonAntielectron";
                    std::string leptonic5 = interacting + "MuonAntitau";
                    if (interactionChannel == leptonic1 || interactionChannel == leptonic2 || interactionChannel == leptonic3 || interactionChannel == leptonic4 || interactionChannel == leptonic5) {
                        this->tabEnergy.erase(this->tabEnergy.begin() + i);
                        this->tabRate.erase(this->tabRate.begin() + i);
                        this->tabE.erase(this->tabE.begin() + i);
                        this->tabs.erase(this->tabs.begin() + i);
                        this->tabCDF.erase(this->tabCDF.begin() + i);
                        this->tabProductsID.erase(this->tabProductsID.begin() + i);
                        
                        sortDictionaryIndexes(i);
                    }
                }
            } else {
                
                int i = el.first;
                
                this->tabEnergy.erase(this->tabEnergy.begin() + i);
                this->tabRate.erase(this->tabRate.begin() + i);
                this->tabE.erase(this->tabE.begin() + i);
                this->tabs.erase(this->tabs.begin() + i);
                this->tabCDF.erase(this->tabCDF.begin() + i);
                this->tabProductsID.erase(this->tabProductsID.begin() + i);
                
                sortDictionaryIndexes(i);
            }
        }
    }
}

std::vector<double> ChannelsBundle::fillTableZeros(std::vector<double> table, size_t size) {
    if (table.size() < size) {
        double difference = size - table.size();
        table.insert(table.begin(), difference, 0);
        return table;
    } else if (table.size() == size) {
        return table;
    } else {
        throw std::runtime_error("Too large table size: suspicious!");
    }
}

void ChannelsBundle::computeInteractionProbabilities (std::vector<std::vector<double>> rates) {
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

void ChannelsBundle::getProductsID(std::vector<double> tabEnergy, double E) {
    // take the index i of the closest tabEnergy[i] to E
    
    double distEnergy = 1e6; // to minimize this value
    int iEnergy = -1;
    for (int i; i <= tabEnergy.size(); ++i) {
        double dist = abs(E - tabEnergy[i]);
        if (dist < distEnergy) {
            distEnergy = dist;
            iEnergy = i;
        }
    }
    
    // select the channel and give the IDs
    std::vector<double> channelProbability = this->channelProbability[iEnergy];
    
    // randomly select the interaction channel
    Random &randomChannel = Random::instance();
    double iRand = randomChannel.rand();
    
    double tabCDFLow = 0;
    double tabCDFUpp = 0;
    int selChannel;
    
    // build a "cumulative function" and select the channel
    for (int i = 0; i <= channelProbability.size() - 1; ++i) {
        tabCDFLow = tabCDFLow + channelProbability[i];
        tabCDFUpp = tabCDFUpp + channelProbability[i + 1];
        if (tabCDFLow < iRand <= tabCDFUpp) {
            selChannel = i;
            break;
        }
    }
    
    std::vector<int> IDs = this->tabProductsID[selChannel];
    this->selectedProductsID = IDs;
    // see the correspondence between the id and the selected channel -> assign the products IDs
    // to take into account parity invariance of the processes... !
    // for the elastic scatterings there are initialising values!
    
    // end of the function to get the IDs!
}

double ChannelsBundle::getRate(int ID, int IDBkg, double E) {
    
    removeChannels(ID, IDBkg);
    
    if (this->tabEnergy.size() == 0) {
        throw std::runtime_error("No active channels for Neutrino-Antineutrino interaction!");
        
    } else if (this->tabEnergy.size() == 1) {
        std::vector<double> ones(this->tabEnergy[0].size(), 1);
        this->channelProbability.push_back(ones); // maybe not needed! Or I need to initialise it in .h
        
        // check if in tabulated energy range
        if (E < this->tabEnergy[0].front() or (E > this->tabEnergy[0].back()))
            throw std::runtime_error("Energy out of range!"); // it is a problem for the simulation! I can put a naive value
        
        getProductsID(this->tabEnergy[0], E);
        
        // interaction rate
        double rate = interpolate(E, this->tabEnergy[0], this->tabRate[0]);
    
        return rate;
        
    } else {
        
        size_t maximumTableSize = -1;
        int maximumIndex = -1;
        for (size_t j = 0; j < this->tabEnergy.size(); j++) {
            size_t tableSize = this->tabRate[j].size();
            if (tableSize > maximumTableSize) {
                maximumTableSize = tableSize;
                maximumIndex = j;
            }
        }
        
        std::vector<double> totalRate(maximumTableSize, 0);
        std::vector<std::vector<double>> resizedVectors;
        
        for (size_t j; j < this->tabRate.size(); ++j) {
            std::vector<double> resizedRate = fillTableZeros(this->tabRate[j], maximumTableSize);
            resizedVectors.push_back(resizedRate);
            
            for (int k = 0; k < maximumTableSize; ++k) {
                totalRate[k] = totalRate[k] + resizedRate[k];
            }
        }
        computeInteractionProbabilities(resizedVectors); // make sure it is sync with the interactionDictionary
        
        // check if in tabulated energy range
        if (E < this->tabEnergy[maximumIndex].front() or (E > this->tabEnergy[maximumIndex].back()))
            throw std::runtime_error("Energy out of range!"); // it is a problem for the simulation! I can put a naive value
        
        getProductsID(this->tabEnergy[maximumIndex], E);
        
        // interaction rate
        double rate = interpolate(E, this->tabEnergy[maximumIndex], totalRate);

        return rate;
    }
}

std::vector<int> ChannelsBundle::getSelectedProductsID() {
    return this->selectedProductsID;
}

} // end namespace nupropa

