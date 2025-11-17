#include "nupropa/Channels.h"
#include <crpropa/Referenced.h>

#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <unordered_map>
#include <iostream>

namespace nupropa {

using namespace crpropa;

Channels::Channels() {};

// maybe we need to code the interaction, as Rhorry did
Channels::Channels(std::vector<std::string> interactionChannels, std::vector<int> active) {
    
    if (interactionChannels.size() != active.size())
        throw std::invalid_argument("Different number of channels and activations!");

    setInteractionChannels(interactionChannels);
    setChannelsActive(active);
    
};

// maybe we need to code the interaction, as Rhorry did
Channels::Channels(std::vector<std::string> interactionChannels, std::vector<int> active, std::string interactionFolderPath) {
    
    if (interactionChannels.size() != active.size())
        throw std::invalid_argument("Different number of channels and activations!");
    
    setInteractionChannels(interactionChannels);
    setChannelsActive(active);
    setInteractionFolderPath(interactionFolderPath);
    
};

int Channels::getChannelIndex(std::string interactionChannel) const {
    int indexChannel = -1;
    
    for (int i; i <= this->interactionChannels.size(); i++) {
        if (this->interactionChannels[i] == interactionChannel) {
            indexChannel = i;
        } else {
            continue;
        }
    }
    return indexChannel;
}

void Channels::setInteractionChannels(std::vector<std::string> interactionChannels) {
    this->interactionChannels = interactionChannels;
}

void Channels::activeAll() {
    std::vector<int> activeAll(this->interactionChannels.size(), true);
    this->active = activeAll;
}

void Channels::setChannelsActive(std::vector<int> active) {
    this->active = active;
}

void Channels::setInteractionFolderPath(std::string interactionFolderPath) {
    this->interactionFolderPath = interactionFolderPath;
}

void Channels::setInactiveChannel(std::string interactionChannel) {
    for (int i; i <= this->interactionChannels.size(); i++) {
        if (this->interactionChannels[i] == interactionChannel) {
            this->active[i] = false;
        } else {
            continue;
        }
    }
}

} // end namespace nupropa

