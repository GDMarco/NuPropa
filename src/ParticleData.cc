#include "nupropa/ParticleID.h"
#include <crpropa/Units.h>

#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <unordered_map>

namespace nupropa {

using namespace crpropa;

ParticleData::ParticleData() {
    
    // to be converted in kg
    std::unordered_map<int, double> IDmass = {
        {211, 0.13957 * GeV / c_squared},   // Pion+
        {111, 0.13498 * GeV / c_squared},   // Pi0
        {321, 0.49368 * GeV / c_squared},   // Kaon+
        {2212, 0.93827 * GeV / c_squared},  // Proton
        {11, 0.000511 * GeV / c_squared},   // Electron
        {13, 0.10566 * GeV / c_squared},    // Muon
        {15, 1.77686 * GeV / c_squared}, // Tauon
        {24, 80.385 * GeV / c_squared}, // W boson
        {23, 91.1876 * GeV / c_squared}, // Z boson
        {1, 0.00467 * GeV / c_squared}, // Down Quark
        {2, 0.00216 * GeV / c_squared}, // Up Quark
        {3, 0.0934 * GeV / c_squared}, // Strange Quark
        {4, 1.24 * GeV / c_squared}, // Charm Quark
        {5, 4.18 * GeV / c_squared}, // Bottom Quark
        {6, 171 * GeV / c_squared} // Top Quark
    }
    
    setParticleIDmass(IDmass);
    
}

ParticleData::ParticleData(std::unordered_map<int, double> IDmass) {
    setParticleIDmass(IDmass);
}

void setParticleIDmass(std::unordered_map<int, double> IDmass) {
    this->dictionaryIDmass = IDmass;
}

double getParticleMassFromID(int ID) {

    auto it = dictionaryIDmass.find(ID);
    if (it != dictionaryIDmass.end()) {
        return it->second;
    } else {
        throw std::runtime_error("Error: Particle ID " << ID << " not found in the dictionary.\n");
        return -1;
    }
}

void addNewParticle(int ID, double mass) {
    if (dictionaryIDmass.find(ID) == dictionaryIDmass.end()) {
        dictionaryIDmass[ID] = mass;
    } else {
        std::cerr << "Warning: Particle ID " << ID << " already exists.\n";
    }
}

} // end namespace nupropa

