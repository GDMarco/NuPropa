#ifndef NUPROPA_CHANNELS_H
#define NUPROPA_CHANNELS_H

#include <crpropa/Referenced.h>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <unordered_map>

namespace nupropa {

using namespace crpropa;

class Channels : public Referenced
{
private:
    
    std::vector<std::string> interactionChannels;
    std::vector<std::vector<int>> productsID; // ID1, ID2(, ID3)
    std::vector<int> active;
    std::string interactionFolderPath = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoAntineutrinoInteraction/";
    
public:
// enum C++
    Channels();
    
    Channels(std::vector<std::string> interactionChannels, std::vector<int> active);
    
    Channels(std::vector<std::string> interactionChannels, std::vector<int> active, std::string interactionFolderPath);
    
    std::vector<std::string> getInteractionChannels() const{
        return this->interactionChannels;
    };
    
    std::vector<int> getActiveChannels() const {
        return this->active;
    };
    std::string getInteractionFolderPath() const {
        return this->interactionFolderPath;
    };
    std::vector<std::vector<int>> getProductsID() const {
        return this->productsID;
    };
    int getChannelIndex(std::string interactioChannel) const;
    
    void activeAll();
    
    void setInteractionChannels(std::vector<std::string> interactionChannels);
    void setChannelsActive(std::vector<int> active);
    void setInteractionFolderPath(std::string interactionFolderPath);
    void setInactiveChannel(std::string interactionChannel);
    
};

} // end namespace nupropa

#endif // CHANNELS_H
