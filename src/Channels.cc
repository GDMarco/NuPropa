#include "Channels.h"

#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>

namespace crpropa {

Channels::Channels(std::vector<std::string> interactionChannels, std::vector<bool> active, std::string interactionFolderPath) {
    if (interactionChannels.size() != active.size())
        throw std::invalid_argument("Different number of channels and activations!")
        // maybe a runtime_error, see documentation!
    
    setInteractionChannels(interactionChannels);
    setChannelsActive(active);
    setInteractionFolderPath(interactionFolderPath);
};

std::vector<std::string> Channels::getInteractionChannels() const {
    return this->interactionChannels;
}

std::vector<bool> Channels::getActiveChannels() const {
    return this->active;
}

std::string getInteractionFolderPath() const {
    return this->interactionFolderPath;
}

std::vector<std::vector<int>> getProductsID() const {
    return this->productsID;
}

int Channels::getChannelIndex(std::string interactionChannel) const {
    int indexChannel = -1;
    
    for (int i, i <= this->interactionChannels.size(), i++) {
        if (this->interactionChannels[i] == interactionChannel) {
            indexChannel = i;
        } else {
            continue
        }
    }
    return indexChannel;
}

void setInteractionChannels(std::vector<std::string> interactionChannels) {
    this->interactionChannels = interactionChannels;
}

std::vector<std::string> Channels::readInteractionChannels(std::string interactionFolderPath) {
    
    std::__fs::filesystem::path dir = interactionFolderPath;
    std::vector<std::string> interactionChannels;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {
        interactionChannels.push_back(dir_entry.path().string());
    }
    
    return interactionChannels;
}

std::vector<std::vector<int>> Channels::readProductsID(std::string interactionFolderPath) {
    std::__fs::filesystem::path dir = interactionFolderPath;
    std::vector<std::string> interactionChannels;
    
    for (auto const& dir_entry : std::__fs::filesystem::directory_iterator{dir}) {
        
        std::string filePath = interactionFolderPath + dir + "products.txt";
        
        // open a productsID.txt file and read the ID numbers
        std::ifstream infile(filePath.c_str());
        if (!infile.good)
            throw std::runtime_error("Could not open " + filePath);
        
        std::string line;
        while (std::getline(infile, line)) {
            std::vector<int> vecIDs;
            
            if ((line.size() > 0) & (line[0] != '#'))
                vecIDs.push_back(std::stoi(line));
        }
        infile.close();
        productsID.push_back(vecIDs);
    }
    return productsID;
}

void Channels::setChannelsActive(std::vector<bool> active) {
    this->active = active;
}

void Channels::setInteractionFolderPath(std::string interactionFolderPath) {
    this->interactionFolderPath = interactionFolderPath;
}

void Channels::setInactiveChannel(std::string interactionChannel) {
    for (int i, i <= this->interactionChannels.size(), i++) {
        if (this->interactionChannels[i] == interactionChannel) {
            this->active[i] = false;
        } else {
            continue;
        }
    }
}

} // end namespace crpropa

