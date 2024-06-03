#ifndef CHANNELS_H
#define CHANNELS_H

#include <crpropa/Referenced.h>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

namespace nupropa {

using namespace crpropa;
/// A custom C++ class for interaction channels
class Channels : public Referenced
{
private:
    
    std::vector<std::string> interactionChannels;
    std::vector<std::vector<int>> productsID; // ID1, ID2(, ID3)
    std::vector<bool> active;
    std::string interactionFolderPath;
    
public:
// enum C++
    Channels(std::vector<std::string> interactionChannels, std::vector<std::vector<int>> productsID, std::vector<bool> active, std::string interactionFolderPath);
    
    std::vector<std::string> getInteractionChannels() const{
        return this->interactionChannels;
    };
    
    std::vector<bool> getActiveChannels() const {
        return this->active;
    };
    std::string getInteractionFolderPath() const {
        return this->interactionFolderPath;
    };
    std::vector<std::vector<int>> getProductsID() const {
        return this->productsID;
    };
    int getChannelIndex(std::string interactioChannel) const;
    
    std::vector<std::string> readInteractionChannels(std::string interactionFolder);
    std::vector<std::vector<int>> readProductsID(std::string interactionFolder);
    
    
    void setInteractionChannels(std::vector<std::string> interactionChannels);
    void setChannelsActive(std::vector<bool> active);
    void setInteractionFolderPath(std::string interactionFolderPath);
    void setProductsID(std::vector<std::vector<int>> productsID);
    void setInactiveChannel(std::string interactionChannel);
};

/**
class NeutrinoAntineutrinoChannels : public Channels {
public:
    std::string interactionFolderPath;
    std::vector<std::string> interactions;
    std::vector<std::vector<int>> productsID;
    std::vector<bool> activeAll;

    NeutrinoAntineutrinoChannels();
};

NeutrinoAntineutrinoChannels::NeutrinoAntineutrinoChannels()
    : Channels(interactions, productsID, activeAll, interactionFolderPath) {
    interactionFolderPath = "/Applications/CRPropa/NuGammaInteraction/CRPropa3-data/data/NeutrinoInteractions/NeutrinoAntineutrinoInteraction/";
    interactions = readInteractionChannels(interactionFolderPath);
    productsID = readProductsID(interactionFolderPath);
    activeAll = std::vector<bool>(interactions.size(), true);
}
 */
} // end namespace nupropa

#endif // CHANNELS_H
