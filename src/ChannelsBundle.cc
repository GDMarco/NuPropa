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

#include <iostream>
#include <sstream>
#include <iomanip>

namespace nupropa {

using namespace crpropa;

ChannelsBundle::ChannelsBundle(ref_ptr<Channels> channels, std::string fname) {
    
    std::string path = channels->getInteractionFolderPath();
    std::vector<std::vector<int>> products = channels->getProductsID();
    std::vector<std::string> interactionChannel = channels->getInteractionChannels();
    std::vector<bool> activeChannels = channels->getActiveChannels();
    
    std::unordered_map<int, std::pair<std::string, std::string >> ratesDict;
    int index = 0;
    
    std::vector<std::string> masses = {"m1", "m2", "m3"};
    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50};
    
    for (int i = 0; i <= interactionChannel.size(); ++i) {
        
        if (!activeChannels[i]) {
            
            continue;
            
        } else {
            
            for (const auto& a : masses) {
                for (const auto& z : redshifts) {
                    
                    std::string pathInt = path + interactionChannel[i];
                    
                    std::ostringstream out;
                    out << std::fixed << std::setprecision(1) << z;
                    std::string zDec = out.str();
                    
                    std::string filename = pathInt + "/rate_" + fname + "_" + a + "_z" + zDec + ".txt";
                    loadRateFile(filename);
                    
                    std::string filenameCDF = pathInt + "/cdf_" + fname + "_" + a + "_z" + zDec + ".txt";
                    loadCumulativeRateFile(filenameCDF);
                    
                    // this is the same for each alphaBeta actually, so a way to not be redundant?
                    std::string filenameProdChan = pathInt + "/products_channelId.txt"; // 2 to 2 processes, so (1 row, 3 col)
                    loadProductsChannelId(filenameProdChan);
                    
                    ratesDict[index] = {interactionChannel[i], fname + "_" + a + "_z" + zDec};
                    index = index + 1;
            
                }
            }
        }
        this->ratesDictionary = ratesDict;
    }
}


void ChannelsBundle::loadRateFile(const std::string& filename) {

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
    
    this->tabEnergy.push_back(vecEnergy);
    this->tabRate.push_back(vecRate);
    
}

void ChannelsBundle::loadCumulativeRateFile(const std::string& filename) {
    
    std::ifstream infile(filename.c_str());
    
    if (!infile.good())
        throw std::runtime_error("Could not open CDF file: " + filename);
    
    std::vector<double> vecE;
    std::vector<double> vecs;
    std::vector<std::vector<double>> vecCDF;
    
    // skip header
    while (infile.peek() == '#')
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    
    double c;
    // read s values in first line
    infile >> c; // skip first value
    while (infile.good() and (infile.peek() != '\n')) {
        infile >> c;
        vecs.push_back(pow(10, c) * eV * eV);
    }
    
    // read all following lines: E, cdf values
    while (infile.good()) {
        infile >> c;
        if (!infile)
            break;  // end of file
        vecE.push_back(pow(10, c) * eV);
        std::vector<double> cdf;
        for (int i = 0; i < vecs.size(); i++) {
            infile >> c;
            cdf.push_back(c / Mpc);
        }
        vecCDF.push_back(cdf);
    }
    infile.close();
    
    this->tabE.push_back(vecE);
    this->tabs.push_back(vecs);
    this->tabCDF.push_back(vecCDF);
    
}

void ChannelsBundle::loadProductsChannelId(const std::string& filename) {
    
    std::ifstream infile(filename.c_str());
    
    if (!infile.good())
        throw std::runtime_error("Could not open rate file: " + filename);
    
    std::vector<int> Ids;
    
    while (infile.good()) {
        if (infile.peek() != '#') {
            double a, b, c;
            infile >> a >> b >> c; // 2 to 2 process, only one row (to optimise)
            if (infile) {
                Ids = {a, b, c};
            }
        }
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infile.close();
    
    this->tabProdChanId.push_back(Ids);
}

std::vector<std::string> ChannelsBundle::getAlphasBetas(int ID, int IDbkg) {
    
    std::vector<std::string> alphasBetas;
    
    if (abs(ID) == abs(IDbkg)) {        // case alpha = beta
        
        std::string interacting = "NeutrinoAntineutrino";
        
        for (const auto& [key, val] : this->ratesDictionary) {
            
            // std::cout << "val.first: " << val.first << '\n';
            std::string interactionChannel = val.first;
            
            if (interactionChannel.find(interacting) == 0) {
                
                // NeutrinoAntineutrinoElastic
                std::string elastic = interacting + "Elastic";
                if (interactionChannel == elastic)
                    alphasBetas.push_back(elastic);
                
                if (abs(ID) == 12) {
                    
                    std::string leptonic = interacting + "Electron";
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                    
                } else if (abs(ID) == 14) {
                    
                    std::string leptonic = interacting + "Muon";
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                    
                } else {
                    
                    std::string leptonic = interacting + "Tauon";
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                    
                } else {
                    
                    alphasBetas.push_back(interactionChannel); // all the other resonances in nu_alpha + nux_alpha
                    
                }
            }
        }
        
    } else {    // for alpha != beta
            
        std::string interacting = "NeutrinoiAntineutrinoj";
            
        for (const auto& [key, val] : ratesDictionary) {
                
            // std::cout << "val.first: " << val.first << '\n';
            std::string interactionChannel = val.first;
            
            if (interactionChannel.find(interacting) == 0) {
                
                // NeutrinoiAntineutrinojElastic
                std::string elastic = interacting + "Elastic";
                if (interactionChannel == elastic)
                    alphasBetas.push_back(elastic);
                
                if ((ID == 12 && IDbkg == -14) || (ID == -14 && IDbkg == 12)) {
                        
                    std::string leptonic = interacting + "ElectronAntimuon";
                        
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
            
                        
                } else if ((ID == 12 && IDbkg == -16) || (ID == -16 && IDbkg == 12)) {
                        
                    std::string leptonic3 = interacting + "ElectronAntitau";
                        
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                        
                } else if ((ID == -12 && IDbkg == 14) || (ID == 14 && IDbkg == -12)) {
                        
                    std::string leptonic = interacting + "MuonAntielectron";
                        
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                    
                } else if ((ID == 14 && IDbkg == -16) || (ID == -16 && IDbkg == 14)) {
                        
                    std::string leptonic = interacting + "MuonAntitau";
                        
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                        
                } else if ((ID == 16 && IDbkg == -12) || (ID == -12 && IDbkg == 16)) {
                        
                    std::string leptonic = interacting + "TauAntielectron";
                        
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                        
                } else if ((ID == 16 && IDbkg == -14) || (ID == -14 && IDbkg == 16)) {
                        
                    std::string leptonic = interacting + "TauAntimuon";
                        
                    if (interactionChannel == leptonic)
                        alphasBetas.push_back(leptonic);
                    
                } else {
                    continue;
                }
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

void ChannelsBundle::computeInteractionProbabilities(std::vector<std::vector<double>> rates) {
    
    this->channelProbability.clear();
    
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

void ChannelsBundle::selectIndex(std::vector<double> tabEnergy, double E) {
    
    // at each step needs to be randomly picked
    this->selectedIndex.clear();
    
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

    this->selectedIndex = selChannel;
    // to take into account parity invariance of the processes... !
    // for the elastic scatterings there are initialising values!
    
}

double ChannelsBundle::findClosestRedshift(double z, const std::vector<double> &redshifts) const {
    
    auto it = std::lower_bound(redshifts.begin(), redshifts.end(), z);

    if (it == redshifts.begin()) return *it;
    if (it == redshifts.end()) return redshifts.back();

    double upper = *it;
    double lower = *(it - 1);

    return (std::abs(upper - z) < std::abs(lower - z)) ? upper : lower;
    
}

void ChannelsBundle::selectIndexes(std::string massCombRedshift, int ID, intID bkg) {
    
    // at each step should be updated
    this->selectedIndexes.clear();
    
    std::vector<std::string> alphasBetas = getAlphasBetas(ID, IDbkg);
    std::vector<int> indexes;

    for (const auto& [key, val] : ratesDictionary) {
        if (std::find(alphasBetas.begin(), alphasBetas.end(), val.first) != alphasBetas.end() && val.second == massCombRedshift) {
            indexes.push_back(key);
        }
    }
    
    this->selectedIndexes = indexes;
}

std::vector<std::vector<double>> ChannelsBundle::selectedRates(std::vector<int> indexes) {
    
    std::vector<std::vector<double>> rates;

    for (int idx : indexes) {
        if (idx >= 0 && idx < this->tabRate.size()) {
            rates.push_back(this->tabRate[idx]);
        }
    }
    
    return rates;
}

std::vector<std::vector<double>> ChannelsBundle::selectedEnergies(std::vector<int> indexes) {
    
    std::vector<std::vector<double>> energies;

    for (int idx : energies) {
        if (idx >= 0 && idx < this->tabEnergy.size()) {
            energies.push_back(this->tabEnergy[idx]);
        }
    }
    
    return energies;
}

std::vector<std::vector<double>> ChannelsBundle::selectCDF() {
    
    std::vector<std::vector<std::vector<double>>> cdfs;

    for (int idx : this->selectedIndexes) {
        if (idx >= 0 && idx < this->tabCDF.size()) {
            cdfs.push_back(this->tabCDF[idx]);
        }
    }
    
    std::vector<std::vector<double>> cdf = cdfs[this->selectedIndex];
    return cdf;
}

std::vector<double> ChannelsBundle::selects(std::vector<int> indexes, int selectedIndex) {
    
    std::vector<std::vector<double>> ss;

    for (int idx : this->selectedIndexes) {
        if (idx >= 0 && idx < this->tabs.size()) {
            ss.push_back(this->tabs[idx]);
        }
    }
    
    std::vector<std::vector<double>> s = ss[this->selectedIndex];
    return s;
}

std::vector<double> ChannelsBundle::selectE() {
    
    std::vector<std::vector<double>> Es;

    for (int idx : this->selectedIndexes) {
        if (idx >= 0 && idx < this->tabE.size()) {
            Es.push_back(this->tabE[idx]);
        }
    }
    
    std::vector<std::vector<double>> E = E[this->selectedIndex];
    return E;
}

std::vector<int> ChannelsBundle::selectProdChanId() {
    
    std::vector<std::vector<int>> prodChanId;

    for (int idx : this->selectedIndexes) {
        if (idx >= 0 && idx < this->tabProdChanId.size()) {
            prodChanId.push_back(this->abProdChanId[idx]);
        }
    }
    
    std::vector<std::vector<double>> IDs = prodChanId[this->selectedIndex];
    return IDs;
    
}

double ChannelsBundle::getRate(int ID, int IDbkg, std::string massComb, double z, double E) {
    
    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50}; // same as computed in NuPropa-data
    double zClosest = findClosestRedshift(z, redshifts);
    std::ostringstream out;
    out << std::fixed << std::setprecision(1) << zClosest;
    std::string zDec = out.str();
    
    std::string redshift = "_z" + zDec;
    
    // look for the indexes compatible with massComb + redshift & alphaBeta
    setSelectedIndexes(massComb + redshift, ID, IDbkg);
    
    std::vector<std::vector<double>> rates = selectedRates(this->selectedIndexes);
    std::vector<std::vector<double>> energies = selectedEnergies(this->selectedIndexes);

    this->channelProbability.clear();
    
    if (energies.size() == 0) {
    
        throw std::runtime_error("No active channels for Neutrino-Antineutrino interaction!");
    
    } else if (energies.size() == 1) { // i.e. indexes has only one element
        std::vector<double> ones(energies.size(), 1);
        this->channelProbability.push_back(ones);
        
        if (E < energies.front() or (E > energies.back()))
            return;
        
        // interaction rate
        double rate = interpolate(E, energies, rates);
        return rate;
        
    } else {
        
        size_t maximumTableSize = -1;
        int maximumIndex = -1;
        for (size_t j = 0; j < energies.size(); ++j) {
            size_t tableSize = rates[j].size();
            if (tableSize > maximumTableSize) {
                maximumTableSize = tableSize;
                maximumIndex = j;
            }
        }
        
        std::vector<double> totalRate(maximumTableSize, 0);
        std::vector<std::vector<double>> resizedRates;
        
        for (size_t j; j < rates.size(); ++j) {
            std::vector<double> resizedRate = fillTableZeros(rates[j], maximumTableSize);
            resizedRates.push_back(resizedRate);
            
            for (int k = 0; k < maximumTableSize; ++k) {
                totalRate[k] = totalRate[k] + resizedRate[k];
            }
        }
        
        computeInteractionProbabilities(resizedRates);
        
        if (E < energies[maximumIndex].front() or (E > energies[maximumIndex].back()))
            return;
        
        getProdChanId(energies[maximumIndex], E);
        
        // interaction rate
        double rate = interpolate(E, this->tabEnergy[maximumIndex], totalRate);
        return rate;
    }
}

std::vector<int> ChannelsBundle::getSelectedIndexes() {
    return this->selectedIndexes;
}

int ChannelsBundle::getSelectedIndex() {
    return this->selectedIndex;
}

} // end namespace nupropa

