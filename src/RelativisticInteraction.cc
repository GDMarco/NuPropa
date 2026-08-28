#include "nupropa/RelativisticInteraction.h"
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <limits>
#include <cmath>
#include <iostream>

namespace nupropa {

using namespace crpropa;

RelativisticInteraction::RelativisticInteraction() {};

RelativisticInteraction::RelativisticInteraction(double m1, double m2, double E, double s) {
    //std::cout << "inside relativistic interaction" << std::endl;
    setBetaCom(E, m1, m2, s);
    setGammaCom(E, s);
};

RelativisticInteraction::RelativisticInteraction(double m1, double m2, double E, double s,
                                                 double epsMin, double epsMax) {
    setBetaCom(E, m1, m2, s, epsMin, epsMax);
    setGammaCom(E, s);
};

RelativisticInteraction::RelativisticInteraction(double m1, double E, double s) {
    setBetaPhotonCom(E, m1, s);
    setGammaCom(E, s);
};

RelativisticInteraction::RelativisticInteraction(double m1, double E, double s,
                                                 double epsMin, double epsMax) {
    setBetaPhotonCom(E, m1, s, epsMin, epsMax);
    setGammaCom(E, s);
};

void RelativisticInteraction::setBetaCom(double E, double m1, double m2, double s,
                                         double /*epsMin*/, double /*epsMax*/) {

    Random &random = Random::instance();
    const double eps = 1e-12;
    const int maxTries = 128;
    const double betaLimitTolerance = 1e-3;
    auto clampBeta = [betaLimitTolerance](double &beta) {
        if (!std::isfinite(beta))
            return false;
        if (std::abs(beta) >= 1.0) {
            if (std::abs(beta) > 1.0 + betaLimitTolerance)
                return false;
            beta = std::copysign(std::nextafter(1.0, 0.0), beta);
        }
        return true;
    };
    for (int t = 0; t < maxTries; ++t) {
        double costh = random.randUniform(-1.0 + eps, 1.0 - eps); // avoid boundaries

        double m1sq = m1 * m1 * c_squared * c_squared;
        double m2sq = m2 * m2 * c_squared * c_squared;
        double p1sq = E * E - m1sq;
        if (p1sq < 0)
            continue;
        double p1 = sqrt(p1sq);
        double K = s - m1sq - m2sq;
        double A = E * E - p1sq * costh * costh;
        if (A <= eps)
            continue;

        double disc = K * K - 4 * A * m2sq;
        if (disc < 0)
            continue;
        double e = (K * E + costh * p1 * sqrt(disc)) / (2 * A);
        if (!std::isfinite(e))
            continue;
        double p2sq = e * e - m2 * m2 * c_squared * c_squared;
        if (p2sq < 0)
            p2sq = 0;

        double beta = (p1 + sqrt(p2sq) * costh) / (E + e);
        if (!clampBeta(beta))
            continue;

        this->eps = e;
        this->beta_com = beta;
        return;
    }

    double m1sq = m1 * m1 * c_squared * c_squared;
    double m2sq = m2 * m2 * c_squared * c_squared;
    double p1sq = E * E - m1sq;
    double K = s - m1sq - m2sq;
    if (p1sq >= 0 && K > 0 && E > 0) {
        double p1 = sqrt(p1sq);
        double targetMassEnergy = sqrt(std::max(m2sq, 0.0));
        double e = 0.0;
        double p2 = 0.0;
        double costh = -1.0;

        if (targetMassEnergy <= 0) {
            e = K / (4 * E);
            p2 = e;
        } else {
            double y = K / (2 * E);
            if (y > 0) {
                e = 0.5 * (y + targetMassEnergy * targetMassEnergy / y);
                double p2sq = e * e - m2sq;
                if (p2sq < 0)
                    p2sq = 0;
                p2 = sqrt(p2sq);
                costh = (y < targetMassEnergy) ? 1.0 : -1.0;
            }
        }

        double beta = (p1 + p2 * costh) / (E + e);
        if (std::isfinite(e) && e > 0 && clampBeta(beta)) {
            this->eps = e;
            this->beta_com = beta;
            return;
        }
    }

    throw std::runtime_error("RelativisticInteraction::setBetaCom failed to find a valid boost");
}

void RelativisticInteraction::setBetaPhotonCom(double E, double m1, double s,
                                               double /*epsMin*/, double /*epsMax*/) {
    Random &random = Random::instance();
    const double eps = 1e-12;
    const int maxTries = 128;
    const double betaLimitTolerance = 1e-3;
    for (int t = 0; t < maxTries; ++t) {
        double costh = random.randUniform(-1.0 + eps, 1.0 - eps); // avoid boundary at +1

        double denom = 1 - costh;
        if (denom < eps)
            denom = eps;
        double e = (s - m1 * m1 * c_squared * c_squared) / 2 / E / denom;
        if (!std::isfinite(e))
            continue;
        double p1sq = E * E - m1 * m1 * c_squared * c_squared;
        if (p1sq < 0)
            p1sq = 0;
        double beta = (sqrt(p1sq) + e * costh) / (E + e);
        if (!std::isfinite(beta))
            continue;
        if (std::abs(beta) >= 1.0) {
            if (std::abs(beta) > 1.0 + betaLimitTolerance)
                continue;
            beta = std::copysign(std::nextafter(1.0, 0.0), beta);
        }

        this->eps = e;
        this->beta_com = beta;
        return;
    }

    throw std::runtime_error("RelativisticInteraction::setBetaPhotonCom failed to find a valid boost");
}

void RelativisticInteraction::setGammaCom(double E, double s) {
    if (!(s > 0) || !std::isfinite(s))
        throw std::runtime_error("RelativisticInteraction::setGammaCom received invalid s");
    double gamma = (E + this->eps) / sqrt(s);
    this->gamma_com = gamma;
}

double RelativisticInteraction::computeProductsMomentumCom(double s, double m3, double m4) {

    double arg = s / 4 - (m3 * m3 * c_squared * c_squared + m4 * m4 * c_squared * c_squared) * 0.5 + std::pow((m3 * m3 * c_squared * c_squared - m4 * m4 * c_squared * c_squared), 2) / 4 / s;

    arg = std::max(arg, 0.0);
    return std::sqrt(arg);
}

std::vector<double> RelativisticInteraction::getProductEnergiesLab(double s, double costh13_com, double m3, double m4) {

    double p_com = computeProductsMomentumCom(s, m3, m4);

    double E3 = this->gamma_com * (sqrt(m3 * m3 * c_squared * c_squared + p_com * p_com) + this->beta_com * p_com * costh13_com);
    double E4 = this->gamma_com * (sqrt(m4 * m4 * c_squared * c_squared + p_com * p_com) - this->beta_com * p_com * costh13_com);

    std::vector<double> energies = {E3, E4};

    for (std::size_t i = 0; i < energies.size(); ++i) {
        if (!std::isfinite(energies[i])) {  // catches both NaN and Inf
            std::cout << "WARNING: Non-finite energy detected!\n";
            std::cout << "  index = " << i << "\n";
            std::cout << "  value = " << energies[i] << "\n";
            std::cout << "  p_com = " << p_com << "\n";
            std::cout << "  m3 = " << m3 << ", m4 = " << m4 << "\n";
            std::cout << "  costh13_com = " << costh13_com << "\n";
            std::cout << "  gamma_com = " << this->gamma_com << "\n";
            std::cout << "  beta_com = " << this->beta_com << "\n";
        }
    }

    return energies;
}

} // end namespace nupropa
