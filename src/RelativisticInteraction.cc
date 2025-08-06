#include "nupropa/RelativisticInteraction.h"
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>

namespace nupropa {

using namespace crpropa;

RelativisticInteraction::RelativisticInteraction() {};

RelativisticInteraction::RelativisticInteraction(double m1, double m2, double E, double s) {
    setBetaCom(E, m1, m2, s);
    setGammaCom(E, s);
};

RelativisticInteraction::RelativisticInteraction(double m1, double E, double s) {
    setBetaPhotonCom(E, m1, s);
    setGammaCom(E, s);
};

void RelativisticInteraction::setBetaCom(double E, double m1, double m2, double s) {
    
    Random &random = Random::instance();
    double costh = random.randUniform(-1, 1); // generate a random number between [-1, 1]
        // it is problematic at the boundaries costh * costh == 1
        // I should take it from the random generation in process!
    
    double K = s - m1 * m1 * c_squared * c_squared - m2 * m2 * c_squared * c_squared;
    double Y = K * K / 4 + m2 * m2 * c_squared * c_squared * E * E * costh * costh;
    
    double e = K * E + sqrt(4 * E * E * (1 - costh * costh) * Y) / 4 / E / E / (1 - costh * costh);
    
    double beta = (sqrt(E * E - m1 * m1 * c_squared * c_squared) + sqrt(e * e - m2 * m2 * c_squared * c_squared) * costh) / (E + e);
    
    // check eps >= m2 * c_squared
    this->eps = e;
    this->beta_com = beta;
    
}

void RelativisticInteraction::setBetaPhotonCom(double E, double m1, double s) {
    Random &random = Random::instance();
    double costh = random.randUniform(-1, 1); // generate a random number between [-1, 1]
    // it should be taken from process!
    
    double e = (s - m1 * m1 * c_squared * c_squared) / 2 / E / (1 - costh);
    
    double beta = (sqrt(E * E - m1 * m1 * c_squared * c_squared) + e * costh) / (E + e);
    
    this->eps = e;
    this->beta_com = beta;
}

void RelativisticInteraction::setGammaCom(double E, double s) {
    double gamma = (E + this->eps) / sqrt(s);
    this->gamma_com = gamma;
}

double RelativisticInteraction::computeProductsMomentumCom(double s, double m3, double m4) {
    
    double p_com = sqrt(s / 4 - (m3 * m3 * c_squared * c_squared + m4 * m4 * c_squared * c_squared) * 0.5 + std::pow((m3 * m3 * c_squared * c_squared - m4 * m4 * c_squared * c_squared), 2) / 4 / s);
    
    return p_com;
}

std::vector<double> RelativisticInteraction::getProductEnergiesLab(double s, double costh13_com, double m3, double m4) {
    
    double p_com = computeProductsMomentumCom(s, m3, m4);
    
    double E3 = this->gamma_com * (sqrt(m3 * m3 * c_squared * c_squared + p_com * p_com) + this->beta_com * p_com * costh13_com);
    double E4 = this->gamma_com * (sqrt(m4 * m4 * c_squared * c_squared + p_com * p_com) - this->beta_com * p_com * costh13_com);
    
    std::vector<double> energies = {E3, E4};
    
    return energies;
}

} // end namespace nupropa

