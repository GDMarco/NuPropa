#include "nupropa/Channels.h"
#include <crpropa/Referenced.h>

#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <unordered_map>

namespace nupropa {

using namespace crpropa;

Channels::Channels() {};

// maybe we need to code the interaction, as Rhorry did
Channels::Channels(std::vector<std::string> interactionChannels, std::vector<std::vector<int>> productsID, std::vector<bool> active, std::string interactionFolderPath) {
    
    if (interactionChannels.size() != active.size())
        throw std::invalid_argument("Different number of channels and activations!");
        // maybe a runtime_error, see documentation!
    
    setInteractionChannels(interactionChannels);
    setProductsID(productsID);
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

void Channels::loadInteractionChannels(std::string interactionFolderPath) {
    
    std::__fs::filesystem::path dir = interactionFolderPath;
    std::vector<std::string> interactionChannels;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {
        interactionChannels.push_back(dir_entry.path().string());
    }
    
    this->interactionChannels = interactionChannels;
}

void Channels::loadProductsID(std::string interactionFolderPath) {
    
    std::__fs::filesystem::path dir = interactionFolderPath;
    std::vector<std::string> interactionChannels;
    std::vector<std::vector<int>> productsID;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {
        
        std::string filePath = dir_entry.path().string() + "/products.txt"; // to check if correct!
        
        // open a productsID.txt file and read the ID numbers
        std::ifstream infile(filePath.c_str());
        if (!infile.good())
            throw std::runtime_error("Could not open " + filePath);
        
        std::string line;
        std::vector<int> vecIDs;
        while (std::getline(infile, line)) {
            if ((line.size() > 0) & (line[0] != '#'))
                vecIDs.push_back(std::stoi(line)); // the first one should be the channel ID, the products ID!
                    // to change in the products.txt file!
            
        }
        infile.close();
        productsID.push_back(vecIDs);
    }
    this->productsID = productsID;
}

void Channels::activeAll() {
    std::vector<bool> activeAll(this->interactionChannels.size(), true);
    this->active = activeAll;
}

void Channels::setChannelsActive(std::vector<bool> active) {
    this->active = active;
}

void Channels::setProductsID(std::vector<std::vector<int>> productsID) {
    this->productsID = productsID;
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

