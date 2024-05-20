#include <crpropa/Referenced.h>

#include <vector>
#include <string>

namespace crpropa {
/// A custom C++ class for interaction channels
class Channels : public Referenced
{
private:
    
    std::vector<std::string> interactionChannels;
    std::vector<bool> active;
    std::string interactionFolderPath;
    
    
public:
// enum C++
    Channels(std::vector<std::string> interactionChannels, std::vector<bool> active, std::string interactionFolderPath);
    
    std::vector<std::string> getInteractionChannels() const;
    std::vector<bool> getActiveChannels() const;
    int getChannelIndex(std::string interactionChannel) const;
    std::string getInteractionFolderPath() const;
    
    std::vector<std::string> readInteractionChannels(std::string interactionFolder);
    
    void setInteractionChannels(std::vector<std::string> interactionChannels);
    void setChannelsActive(std::vector<bool> active);
    void setInteractionFolderPath(std::string interactionFolderPath);
    
    void setInactiveInteraction(std::string interactionChannel);
};


class NeutrinoAntineutrinoChannels: public Channels {
public:
    std::string interactionFolderPath = "/Applications/CRPropa/NuGammaInteraction/CRPropa3-data/data/NeutrinoAntineutrinoInteraction/"
    std::vector<std::string> interactions = readInteractionChannels(interactionFolderPath);
    std::vector<bool> activeAll(interactions.size(), true);
    
    NeutrinoAntineutrinoChannels() : Channels(interactions, activeAll, interactionFolderPath);
}


} // end namespace crpropa
