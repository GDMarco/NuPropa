#ifndef NUPROPA_NEUTRINOMIXING_H
#define NUPROPA_NEUTRINOMIXING_H

#include <crpropa/Module.h>
#include <crpropa/Referenced.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>


#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace nupropa {

using namespace crpropa;
/// A custom C++ module for flavour-mass conversion and viceversa
class NeutrinoMixing : public Referenced
{
private:
    
    double theta13, theta23, theta12, delta; // parameters for the U matrix
    
    double mass1 = 0;
    double mass2 =  8.3e-3;
    double mass2 =  50e-3; // eV, maybe needed initialisation
    
    Eigen::Matrix3cd UpmnsMatrix;
    Eigen::Matrix3cd flavourMassProbabilities;
    
    std::array<double, 3> massValues = {mass1, mass2, mass3};
    
    double ReU4[3][3][3];
    double ImU4[3][3][3];
    
public:
    
    NeutrinoMixing();
    NeutrinoMixing(double theta13, double theta23, double theta12, double delta);
    
    void setTheta(double theta13, double theta23, double theta12) {
        this->theta13 = theta13;
        this->theta23 = theta23;
        this->theta12 = theta12;
    }
    void setDelta(double delta) {
        this->delta = delta;
    }
    void setMasses(double mass1, double mass2, double mass3) {
        this->massValues = {mass1, mass2, mass3};
        // check for consistency with experimental bounds on the sum
    }
    
    Eigen::Matrix3cd buildUpmnsMatrix();
    
    Eigen::Matrix3cd getUpmnsMatrix() {
        return this->UpmnsMatrix;
    }
    
    Eigen::Matrix3cd getFlavourMassProbabilities() {
        return this->flavourMassProbabilities;
    }
    
    int IdToFlavorIndex(int ID);
    int flavorIndexToId(int index);
    double massIndexToMass(int index); // mass returned in eV
    int massToMassIndex(double mass);
    
    double fromFlavourToMass(int ID);
    int fromMassToFlavour(double mass);
    
    void precomputeOscillationTerms();
    int oscillateFlavour(int ID, double E, double L);
    
};

} // end namespace nupropa

#endif
