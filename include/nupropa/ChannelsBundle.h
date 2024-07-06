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
/// A custom C++ class for interaction channels
class ChannelsBundle : public Referenced
{
private:
    
    std::vector<std::vector<double>> tabEnergy; // 4 columns table depending on the neutrino flavour and number, in initRate they should be built properly
    std::vector<std::vector<double>> tabRate;
    
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    
    std::vector<std::vector<int>> tabProductsID;
    std::vector<int> selectedProductsID;
    
    std::unordered_map<int, std::string> interactionDictionary;
    std::vector<std::vector<double>> channelProbability; // to be sync with interactionDictionary
    
public:
 
    ChannelsBundle();
    ChannelsBundle(ref_ptr<Channels> channels, std::string fname);
    
    double getRate(int ID, int IDBkg, double E);
    void removeChannels(int ID, int IDBkg);
    void sortDictionaryIndexes(int indexErased);
    std::vector<double> fillTableZeros(std::vector<double> table, size_t size);
    void computeInteractionProbabilities(std::vector<std::vector<double>> rates);
    void getProductsID(std::vector<double> tabEnergy, double E);
    std::vector<int> getSelectedProductsID();
    
};

} // end namespace nupropa

#endif // NUPROPA_CHANNELSBUNDLE_H
