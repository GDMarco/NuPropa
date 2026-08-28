#include "nupropa/Hadronisations.h"

#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Vector3.h"
#include "crpropa/Candidate.h"

#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Pythia8/Pythia.h"

namespace Pythia8 {

class PythiaHadronisation {
private:

    std::vector<std::vector<double>> secondaries;

public:
    PythiaHadronisation(std::vector<int> Ids, std::vector<double> Es, crpropa::Vector3d dir, std::vector<int> IdDecaying) {
        std::string xmlDir = "./";
        bool printBanner = false;
        Pythia pythia(xmlDir, printBanner); // pythia instance

        int seed = crpropa::Random::instance().randInt(900000000);
        pythia.readString("Random:setSeed = on");
        pythia.readString("Random:seed = " + std::to_string(seed));

        // Hadronise a manually supplied external partonic state. Pythia's full
        // event check expects a beam/process record, so validate the resulting
        // final state explicitly below instead.
        pythia.readString("ProcessLevel:all = off");
        pythia.readString("PartonLevel:all = off");
        pythia.readString("HadronLevel:all = on");
        pythia.readString("HadronLevel:Hadronize = on");
        pythia.readString("Check:event = off");
        pythia.readString("Print:quiet = on");
        pythia.readString("Init:showAllSettings = off");
        pythia.readString("Init:showChangedSettings = off");
        pythia.readString("Init:showAllParticleData = off");
        pythia.readString("Init:showChangedParticleData = off");
        pythia.readString("Next:numberShowInfo = 0");
        pythia.readString("Next:numberShowProcess = 0");
        pythia.readString("Next:numberShowEvent = 0");

        for (std::size_t i = 0; i < IdDecaying.size(); i++)
            pythia.readString(std::to_string(IdDecaying[i]) + ":mayDecay = on"); // enable particle decay

        pythia.readString("Print:quiet = on");
        pythia.readString("Init:showAllSettings = off");
        pythia.readString("Init:showAllParticleData = off");

        if (!pythia.init())
            throw std::runtime_error("PYTHIA initialization failed");
        generateSecondaries(Ids, Es, dir, pythia);

    };

    double getMomentum(double E, double mass) {
        const double thresholdTolerance = 1e-4;
        if (E < mass && (mass - E) / std::max(mass, 1.0) <= thresholdTolerance)
            return 0.0;
        if (E < mass)
            throw std::runtime_error("Energy must be larger than mass");

        return sqrt(std::max(0.0, E * E - mass * mass));
    };

    void generateSecondaries(std::vector<int> Id, std::vector<double> E, crpropa::Vector3d dir, Pythia& pythia) {
        if (Id.empty() || Id.size() != E.size())
            throw std::runtime_error("Invalid PYTHIA input particle lists");

        this->secondaries.clear();

        double mass = 0;
        double p = 0;
        Event& event = pythia.event;
        event.reset();

        double sumPx = 0.0;
        double sumPy = 0.0;
        double sumPz = 0.0;
        double sumE = 0.0;

        for (std::size_t i = 0; i < Id.size(); ++i) {
            if (!std::isfinite(E[i]) || E[i] <= 0)
                throw std::runtime_error("Invalid PYTHIA input energy");
            mass = pythia.particleData.m0(Id[i]);
            p = getMomentum(E[i], mass);
            crpropa::Vector3d d = dir;
            if (Id.size() == 2 && i == 1)
                d = crpropa::Vector3d(-dir.x, -dir.y, -dir.z);
            double px = p * d.x;
            double py = p * d.y;
            double pz = p * d.z;
            sumPx += px;
            sumPy += py;
            sumPz += pz;
            sumE += E[i];
            int col = 0;
            int acol = 0;
            if (std::abs(Id[i]) <= 6) {
                if (Id[i] > 0) {
                    col = 101;
                    acol = 0;
                } else {
                    col = 0;
                    acol = 101;
                }
            }
            event.append(Id[i], 23, col, acol, px, py, pz, E[i], mass);
        }

        if (!pythia.forceHadronLevel(false))
            throw std::runtime_error("PYTHIA forced hadronisation failed");

        std::vector<std::vector<double>> generatedSecondaries;
        double finalEnergy = 0.0;
        double finalPx = 0.0;
        double finalPy = 0.0;
        double finalPz = 0.0;
        for (std::size_t i = 0; i < event.size(); ++i) {
            if (event[i].isFinal()) {
                if (!std::isfinite(event[i].e()) || event[i].e() < 0 ||
                    !std::isfinite(event[i].px()) || !std::isfinite(event[i].py()) ||
                    !std::isfinite(event[i].pz())) {
                    throw std::runtime_error("PYTHIA produced an invalid final-state four-momentum");
                }
                std::vector<double> prop = { (double) event[i].id(), event[i].px(), event[i].py(), event[i].pz(), event[i].e() };
                generatedSecondaries.push_back(prop);
                finalEnergy += event[i].e();
                finalPx += event[i].px();
                finalPy += event[i].py();
                finalPz += event[i].pz();
            }
        }

        const double energyTolerance = 1e-2;
        const double energyScale = std::max(std::abs(sumE), 1.0);
        const double relativeMismatch = std::abs(finalEnergy - sumE) / energyScale;
        const double momentumMismatch = std::sqrt(
            std::pow(finalPx - sumPx, 2) +
            std::pow(finalPy - sumPy, 2) +
            std::pow(finalPz - sumPz, 2)) / energyScale;
        if (!std::isfinite(finalEnergy) || !std::isfinite(momentumMismatch) ||
            relativeMismatch > energyTolerance || momentumMismatch > energyTolerance) {
            std::ostringstream message;
            message << "PYTHIA final-state four-momentum mismatch: input=" << sumE
                    << " GeV output=" << finalEnergy
                    << " GeV energyMismatch=" << relativeMismatch
                    << " momentumMismatch=" << momentumMismatch;
            throw std::runtime_error(message.str());
        }

        this->secondaries = std::move(generatedSecondaries);
    };

    std::vector<std::vector<double>> getSecondaries() {
        return this->secondaries;
    };
};
} // end Pythia8 namespace

namespace nupropa {

using namespace crpropa;

Hadronisations::Hadronisations(bool haveOtherSecondaries, bool haveNeutrinos, bool angularCorrection) { // double thinning,

    // setThinning(thinning);
    setHaveOtherSecondaries(haveOtherSecondaries);
    setHaveNeutrinos(haveNeutrinos);
    setAngularCorrection(angularCorrection);
    setDescription("Hadronisation from PYTHIA");

}

void Hadronisations::setHaveOtherSecondaries(bool otherSecondaries) {
    this->haveOtherSecondaries = otherSecondaries;
}

void Hadronisations::setHaveNeutrinos(bool neutrinos) {
    this->haveNeutrinos = neutrinos;
}

void Hadronisations::setAngularCorrection(bool correction) {
    this->angularCorrection = correction;
}

/**
 void Hadronisations::setThinning(double thinning) {
 this->thinning = thinning;
 }
 */

void Hadronisations::performHadronisation(Candidate *candidate, std::vector<int> Ids, std::vector<double> Es, Vector3d pos, double w) const {
    if (Ids.empty() || Ids.size() != Es.size()) {
        std::cout << "Hadronisation rejected: invalid input particle lists" << std::endl;
        return;
    }

    for (double energy : Es) {
        if (!std::isfinite(energy) || energy <= 0) {
            std::cout << "Hadronisation rejected: invalid input energy" << std::endl;
            return;
        }
    }

    std::vector<int> IdDecaying = {13, 15, 23, 24, 211, 111, 130, 310, 311, 321};
    Vector3d dir = candidate->current.getDirection();
    if (!std::isfinite(dir.x) || !std::isfinite(dir.y) || !std::isfinite(dir.z) ||
        dir.getR() == 0) {
        std::cout << "Hadronisation rejected: invalid input direction" << std::endl;
        return;
    }
    dir = dir / dir.getR();

    std::vector<double> energiesGeV;
    energiesGeV.reserve(Es.size());
    for (double e : Es) {
        energiesGeV.push_back(e / GeV);
    }

    std::vector<std::vector<double>> secondaries;

    const int maxPythiaAttempts = 3;
    bool validHadronisation = false;
    std::string lastPythiaError;
    for (int attempt = 0; attempt < maxPythiaAttempts && !validHadronisation; ++attempt) {
#pragma omp critical(PYTHIAevent)
        {
            try {
                Pythia8::PythiaHadronisation pythiaHadronisation(Ids, energiesGeV, dir, IdDecaying);
                secondaries = pythiaHadronisation.getSecondaries();
                validHadronisation = true;
            } catch (const std::exception& ex) {
                lastPythiaError = ex.what();
            }
        }
    }

    if (!validHadronisation) {
        double inputEnergy = 0.0;
        for (double energy : Es)
            inputEnergy += energy;
        std::cout << "Hadronisation rejected after " << maxPythiaAttempts
                  << " attempts: " << lastPythiaError
                  << " input=" << inputEnergy / EeV << " EeV"
                  << std::endl;
        return;
    }

    std::string tag = getHadronisationTag();

    for (const std::vector<double>& row : secondaries) {
        int id = int(row[0]);
        bool isNeutrino = (std::abs(id) == 12 || std::abs(id) == 14 || std::abs(id) == 16);
        if ((isNeutrino && !haveNeutrinos) || (!isNeutrino && !haveOtherSecondaries))
            continue;

        // Preserve CRPropa's standard SN1/ID1/E1 parent-state bookkeeping.
        candidate->addSecondary(id, row[4] * GeV, pos, w, tag);

        if (angularCorrection) {
            Vector3d secondaryDirection(row[1], row[2], row[3]);
            if (secondaryDirection.getR() > 0) {
                candidate->secondaries.back()->current.setDirection(
                    secondaryDirection / secondaryDirection.getR());
            }
        }
    }
}

void Hadronisations::setHadronisationTag(std::string tag) const {
    this->hadronisationTag = tag;
}

std::string Hadronisations::getHadronisationTag() const {
    return this->hadronisationTag;
}

void Hadronisations::setPendingHadronisation(std::vector<int> Ids, std::vector<double> Es, Vector3d pos, double w) const {
    this->pendingIds = std::move(Ids);
    this->pendingEs = std::move(Es);
    this->pendingPos = pos;
    this->pendingW = w;
    this->hasPending = true;
}

void Hadronisations::process(Candidate *candidate) const {
    if (!hasPending)
        return;
    performHadronisation(candidate, pendingIds, pendingEs, pendingPos, pendingW);
    hasPending = false;
}

} // end namespace nupropa
