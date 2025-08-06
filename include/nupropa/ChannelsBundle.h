#ifndef NUPROPA_CHANNELSBUNDLE_H
#define NUPROPA_CHANNELSBUNDLE_H

#include "nupropa/Channels.h"

#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

namespace nupropa {

using namespace crpropa;

class ChannelsBundle : public Referenced
{
private:
    
    std::vector<std::vector<double>> tabEnergy; // 4 columns table depending on the neutrino flavour and number, in initRate they should be built properly
    std::vector<std::vector<double>> tabRate;
    
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    std::vector<std::vector<int>> tabProdChanId;
    
    std::vector<int> selectedIndexes; // 1st step of selecting the rate
    int selectedIndex; // 2nd step to take the channel for the perform interaction
    
    std::unordered_map<int, std::pair<std::string, std::string>> ratesDictionary;
    std::vector<std::vector<double>> channelProbability;
    
public:
 
    ChannelsBundle(ref_ptr<Channels> channels, std::string fname);
    
    void loadRateFile(const std::string& filename);
    void loadCumulativeRateFile(const std::string& filename);
    void loadProductsChannelId(const std::string& filename);
    
    double findClosestRedshift(double z, const std::vector<double> &redshifts) const;
    void selectIndexes(std::string massCombRedshift, int ID, int IDbkg);
    std::vector<std::string> getAlphasBetas(int ID, int IDbkg) ;
    std::vector<int> selectProdChanId();
    
    std::vector<std::vector<double>> selectedRates(const std::vector<int>& indexes);
    std::vector<std::vector<double>> selectedEnergies(const std::vector<int>& indexes);
    
    double getRate(int ID, int IDBkg, std::string massComb, double z, double E);
    
    std::vector<double> fillTableZeros(std::vector<double> table, size_t size);
    void computeInteractionProbabilities(std::vector<std::vector<double>> rates);
    void selectIndex(std::vector<double> tabEnergy, double E);
    
    std::vector<std::vector<double>> selectCDF();
    std::vector<double> selects();
    std::vector<double> selectE();

    std::vector<int> getSelectedIndexes();
    int getSelectedIndex();
    
};

} // end namespace nupropa

#endif // NUPROPA_CHANNELSBUNDLE_H
