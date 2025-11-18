#include "nupropa/NeutrinoMixing.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/ParticleData.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <string>
#include <complex>
#include <filesystem>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <cmath>

namespace nupropa {

using namespace crpropa;

NeutrinoMixing::NeutrinoMixing() {};

NeutrinoMixing::NeutrinoMixing(double theta13, double theta23, double theta12, double delta) {
    
    setTheta(theta13, theta23, theta12);
    setDelta(delta);
    
    buildUpmnsMatrix();
    precomputeOscillationTerms();
    
}

void NeutrinoMixing::buildUpmnsMatrix() {
        
    double s12 = std::sin(theta12), c12 = std::cos(theta12);
    double s13 = std::sin(theta13), c13 = std::cos(theta13);
    double s23 = std::sin(theta23), c23 = std::cos(theta23);
    std::complex<double> expNegDelta = std::exp(std::complex<double>(0, -delta));
    std::complex<double> expDelta = std::exp(std::complex<double>(0, delta));

    Eigen::Matrix3cd U;
    Eigen::Matrix3d Usquared;
    
    U(0, 0) = c12 * c13;
    U(0, 1) = s12 * c13;
    U(0, 2) = s13 * expNegDelta;

    U(1, 0) = -s12 * c23 - c12 * s23 * s13 * expDelta;
    U(1, 1) = c12 * c23 - s12 * s23 * s13 * expDelta;
    U(1, 2) = s23 * c13;

    U(2, 0) =  s12 * s23 - c12 * c23 * s13 * expDelta;
    U(2, 1) = -c12 * s23 - s12 * c23 * s13 * expDelta;
    U(2, 2) = c23 * c13;
    
    for (int alpha = 0; alpha < 3; ++alpha) {
        for (int i = 0; i < 3; ++i) {
            Usquared(alpha, i) = std::norm(U(alpha, i)); // |U_{αi}|²
        }
    }
    
    // check and validation for unitarity?
    this->UpmnsMatrix = U;
    this->flavourMassProbabilities = Usquared;
    
    precomputeOscillationTerms();
    
}
 
int NeutrinoMixing::flavourIndexToId(int index) {
    static const int pdgTable[3] = {12, 14, 16};
    if (index < 0 || index > 2)
        return 0;

    int pdg = pdgTable[index];
    return pdg;
}

int NeutrinoMixing::IdToFlavourIndex(int ID) {
    static const std::unordered_map<int, int> pdgToFlavourIndex = {
        {12, 0}, {-12, 0},  // ν_e, ν̄_e
        {14, 1}, {-14, 1},  // ν_μ, ν̄_μ
        {16, 2}, {-16, 2}   // ν_τ, ν̄_τ
    };

    auto it = pdgToFlavourIndex.find(ID);
    return (it != pdgToFlavourIndex.end()) ? it->second : -1;
}

double NeutrinoMixing::massIndexToMass(int index) {
    // masses in eV
    static const double massTable[3] = {
        0.0,        // m1
        8.3e-3,     // m2
        5.0e-2      // m3
    };

    if (index < 0 || index > 2)
        return -1.0; // invalid

    return massTable[index];
}

int NeutrinoMixing::massToIndexMass(double mass) {
    static const double epsilon = 1e-4; // tolerance for float comparison
    
    for (int i = 0; i < 3; ++i) {
        if (std::abs(mass - this->massValues[i]) < epsilon)
            return i;
    }
    return -1; // not found
}

double NeutrinoMixing::fromFlavourToMass(int ID) {
    
    int alpha = IdToFlavourIndex(ID);
    if (alpha < 0) return -1;
    
    double r = crpropa::Random::instance().rand();
    double cumulative = 0.0;
    int indexMass = -1;
    
    for (int i = 0; i < 3; ++i) {
        cumulative = cumulative + flavourMassProbabilities(alpha, i);
        if (r < cumulative) {
            indexMass = i;
            return massIndexToMass(indexMass); // in eV
        }
    }
    return -2;
}

int NeutrinoMixing::fromMassToFlavour(double mass) {
    
    int massIndex = massToIndexMass(mass);
    
    if (massIndex < 0)
        return -1;

    double r = crpropa::Random::instance().rand();
    double cumulative = 0;
    int index;
    
    for (int alpha = 0; alpha < 3; ++alpha) {
        cumulative = cumulative + flavourMassProbabilities(alpha, massIndex); // |U_{α i}|²
        if (r < cumulative) {
            index = alpha;
            return flavourIndexToId(index);
        }
    }
    return -2;
}

void NeutrinoMixing::precomputeOscillationTerms() {
    for (int alpha = 0; alpha < 3; ++alpha) {
        for (int beta = 0; beta < 3; ++beta) {
            // (i=1, j=0)
            std::complex<double> term10 = UpmnsMatrix(alpha, 1) * std::conj(UpmnsMatrix(beta, 1)) * std::conj(UpmnsMatrix(alpha, 0)) * UpmnsMatrix(beta, 0);
            ReU4[alpha][beta][0] = term10.real();
            ImU4[alpha][beta][0] = term10.imag();
            
            // (i=2, j=0)
            std::complex<double> term20 = UpmnsMatrix(alpha, 2) * std::conj(UpmnsMatrix(beta, 2)) * std::conj(UpmnsMatrix(alpha, 0)) * UpmnsMatrix(beta, 0);
            
            ReU4[alpha][beta][1] = term20.real();
            ImU4[alpha][beta][1] = term20.imag();
            
            // (i=2, j=1)
            std::complex<double> term21 = UpmnsMatrix(alpha, 2) * std::conj(UpmnsMatrix(beta, 2)) * std::conj(UpmnsMatrix(alpha, 1)) * UpmnsMatrix(beta, 1);
            ReU4[alpha][beta][2] = term21.real();
            ImU4[alpha][beta][2] = term21.imag();
        }
    }
}

int NeutrinoMixing::oscillateFlavour(int ID, double E, double L) {
    
    int alpha = IdToFlavourIndex(ID);
    int sign = (ID < 0) ? -1 : +1;
    
    double conversionPhase = eV * eV * eV * eV / h_planck * 2 * M_PI / c_light / 2;
    
    double phase10 = (massValues[1] * massValues[1] - massValues[0] * massValues[0]) * L / E;
    double phase20 = (massValues[2]  * massValues[2] - massValues[0] * massValues[0]) * L / E * conversionPhase;
    double phase21 = (massValues[2] * massValues[2] - massValues[1] * massValues[1]) * L / E * conversionPhase;

    double P[3]; // oscillation probabilities α → β
    for (int beta = 0; beta < 3; ++beta) {
        double prob = (alpha == beta) ? 1.0 : 0.0;
        
        // real (CP-even) terms
        prob -= 4 * ReU4[alpha][beta][0] * std::pow(std::sin(phase10 / 2), 2);
        prob -= 4 * ReU4[alpha][beta][1] * std::pow(std::sin(phase20 / 2), 2);
        prob -= 4 * ReU4[alpha][beta][2] * std::pow(std::sin(phase21 / 2), 2);
        
        // imaginary terms
        prob += sign * 2 * ImU4[alpha][beta][0] * std::sin(phase10);
        prob += sign * 2 * ImU4[alpha][beta][1] * std::sin(phase20);
        prob += sign * 2 * ImU4[alpha][beta][2] * std::sin(phase21);
        
        P[beta] = clip(prob, 0.0, 1.0);
    }
    
    // normalize (in case of small numerical drift) TO CHECK!
    double total = P[0] + P[1] + P[2];
    for (int i = 0; i < 3; ++i)
        P[i] /= total;
    
    // randomly select oscillated flavour index
    double r = crpropa::Random::instance().rand();
    double cumulative = 0;
    int IdOscillated;
    for (int beta = 0; beta < 3; ++beta) {
        cumulative += P[beta];
        
        if (r < cumulative) {
            IdOscillated = sign * flavourIndexToId(beta);
            return IdOscillated;
        }
    }
}

} // end namespace nupropa

