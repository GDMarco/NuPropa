#include "nupropa/NeutrinoNeutrinoInteraction.h"
#include "nupropa/NeutrinoField.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/NeutrinoMixing.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>


#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace nupropa {

using namespace crpropa;

NeutrinoNeutrinoInteraction::NeutrinoNeutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, ref_ptr<NeutrinoMixing> neutrinoMixing, bool haveSecondaries, double limit) : Module(), secondariesDistribution("dsigdcosth") {

    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    setNeutrinoField(neutrinoField);
    setNeutrinoMixing(neutrinoMixing);

}

void NeutrinoNeutrinoInteraction::setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField) {
    this->neutrinoField = neutrinoField;
    this->neutrinoFieldMass = neutrinoField->getMass();
    std::string fname = neutrinoField->getFieldName();
    setDescription("NeutrinoNeutrinoInteraction::Module" + fname);
    setInteractionTag("NuNuInt");

    std::string pathNuNu = "/sdf/home/g/gaetano/CRPropa/CRPropa3-data/dataOff/NeutrinoInteractions/NeutrinoNeutrinoInteraction/NeutrinoNeutrinoElastic/";
    // getDataPath("NeutrinoNeutrinoInteraction/NeutrinoNeutrinoElastic/"
    std::string pathNuiNuj = "/sdf/home/g/gaetano/CRPropa/CRPropa3-data/dataOff/NeutrinoInteractions/NeutrinoNeutrinoInteraction/NeutrinoiNeutrinojElastic/";
    // getDataPath("NeutrinoNeutrinoInteraction/NeutrinoiNeutrinojElastic/"

    initRate(pathNuNu + "rate_" + fname, pathNuiNuj + "rate_" + fname);
    initCumulativeRate(pathNuNu + "cdf_" + fname, pathNuiNuj + "cdf_" + fname);

    buildSecondariesDistributionClass();


}

void NeutrinoNeutrinoInteraction::setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    this->neutrinoMixing = neutrinoMixing;
}

void NeutrinoNeutrinoInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoNeutrinoInteraction::setLimit(double limit) {
    this->limit = limit;
}

void NeutrinoNeutrinoInteraction::loadRateFile(const std::string& fileName) {
    std::ifstream infile(fileName.c_str());
    if (!infile.good()) {
        throw std::runtime_error("Could not open file: " + fileName);
    }

    std::vector<double> vecEnergy, vecRate;

    while (infile.good()) {
        if (infile.peek() != '#') {
            double a, b;
            infile >> a >> b;
            if (infile) {
                vecEnergy.push_back(std::pow(10, a) * eV);
                vecRate.push_back(b / Mpc);
            }
        }
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    infile.close();

    this->tabEnergy.push_back(vecEnergy);
    this->tabRate.push_back(vecRate);
}

void NeutrinoNeutrinoInteraction::initRate(std::string pathNuNu, std::string pathNuiNuj) {

    this->tabEnergy.clear();
    this->tabRate.clear();

    std::vector<std::string> masses = {"m1", "m2", "m3"};
    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50}; // same as computed in NuPropa-data

    int i = 0;
    std::unordered_map<std::string, int> ratesDict;

    // 9 combination of masses! m1_m1, m1_m2, m2_m1, ...
    /**
     m1_m1_z0 NuNu 0
     m1_m1_z0 NuiNuj 1
     m1_m2_z0 NuNu 2
    ...
     */
    for (const auto& a : masses) {
        for (const auto& z : redshifts) {

            std::ostringstream out;
            out << std::fixed << std::setprecision(1) << z;
            std::string zDec = out.str();

            loadRateFile(pathNuNu + "_" + a + "_z" + zDec + ".txt");
            ratesDict[this->neutrinoField->getFieldName() + "_" + a + "_z" + zDec + "_NuNu"] = i;
            i = i + 1;

            loadRateFile(pathNuiNuj + "_" + a + "_z" + zDec + ".txt");
            ratesDict[this->neutrinoField->getFieldName() + "_" + a + "_z" + zDec + "_NuiNuj"] = i;
            i = i + 1;
        }
    }
    this->ratesDictionary = ratesDict;
}

void NeutrinoNeutrinoInteraction::loadCumulativeRateFile(const std::string& fileName) {

    std::ifstream infile(fileName.c_str());

    if (!infile.good())
        throw std::runtime_error("NeutrinoNeutrinoInteraction: could not open file" + fileName);

    // skip header
    while (infile.peek() == '#')
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');

    // declare the row vectors
    std::vector<double> vecE;
    std::vector<double> vecs;
    std::vector<std::vector<double>> vecCDF;

    // read s values in the first line
    double a;
    infile >> a; // skip the first value
    while (infile.good() and (infile.peek() != '\n')) {
        infile >> a;
        vecs.push_back(pow(10, a) * eV * eV);
    }

    // read all the following lines: E, cdf values
    while(infile.good()) {
        infile >> a;
        if (!infile)
            break; // end of file

        vecE.push_back(pow(10, a) * eV);

        std::vector<double> cdf;
        for (int i = 0; i < vecs.size(); i++) {
            infile >> a;
            cdf.push_back(a / Mpc);
        }
        vecCDF.push_back(cdf);
    }

    this->tabs.push_back(vecs);
    this->tabE.push_back(vecE);
    this->tabCDF.push_back(vecCDF);

}

void NeutrinoNeutrinoInteraction::initCumulativeRate(std::string pathNuNu, std::string pathNuiNuj) {

    this->tabE.clear();
    this->tabs.clear();
    this->tabCDF.clear();

    std::vector<std::string> masses = {"m1", "m2", "m3"};
    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50}; // same as computed in NuPropa-data

    for (const auto& a : masses) {
        for (const auto& z : redshifts) {

            std::ostringstream out;
            out << std::fixed << std::setprecision(1) << z;
            std::string zDec = out.str();

            loadCumulativeRateFile(pathNuNu + "_" + a + "_z" + zDec + ".txt");
            loadCumulativeRateFile(pathNuiNuj + "_" + a + "_z" + zDec + ".txt");
        }
    }
}

void NeutrinoNeutrinoInteraction::setRelativisticInteraction(double m1, double m2, double E, double s) const {
    this->relInteraction = new RelativisticInteraction(m1, m2, E, s);
}

void NeutrinoNeutrinoInteraction::buildSecondariesDistributionClass() {

    if (!haveSecondaries)
        return;

    std::vector<int> channelsId = {1, 3, 5, 8};

    for (int i = 0; i < channelsId.size(); ++i) {
        this->secondariesDistribution.buildChannel(channelsId[i]);
    }
}

void NeutrinoNeutrinoInteraction::performInteraction(Candidate *candidate, int index, double mass, int IDbkg) const {

    if (!haveSecondaries) {
        candidate->setActive(false);
        return;
    }

    double E = candidate->current.getEnergy();
    int ID = candidate->current.getId();
    double w = 1; // no thinning, TBD

    // here I just need the index found before in computeRate
    std::vector<double> vecE = this->tabE[index];
    std::vector<double> vecs = this->tabs[index];
    std::vector<std::vector<double>> vecCDF = this->tabCDF[index];

    // check if in tabulated energy range
    if (vecE.empty() || vecs.empty() || vecCDF.empty()) {
        std::cout << "NuNu no secondaries: empty tables"
                  << " index=" << index << std::endl;
        return;
    }
    if (E < vecE.front() or (E > vecE.back())) {
        std::cout << "NuNu no secondaries: E outside secondary table"
                  << " E=" << E / EeV
                  << " Emin=" << vecE.front() / EeV
                  << " Emax=" << vecE.back() / EeV
                  << std::endl;
            return;
    }

    // sample the value of s
    Random &random = Random::instance();
    size_t i = closestIndex(E, vecE);  // find closest tabulation point
    if (i >= vecCDF.size()) {
        std::cout << "NuNu no secondaries: CDF index out of range"
                  << " i=" << i << " size=" << vecCDF.size() << std::endl;
        return;
    }
    size_t j = random.randBin(vecCDF[i]);
    if (j >= vecs.size()) {
        std::cout << "NuNu no secondaries: s-bin out of range"
                  << " j=" << j << " size=" << vecs.size() << std::endl;
        return;
    }
    double hi = vecs[j];
    double sThr = std::pow(mass + this->neutrinoFieldMass, 2);
    double lo = (j == 0) ? sThr : std::max(sThr, vecs[j - 1]); // first s-tabulation point below min(s_kin)
    if (hi < lo)
        hi = lo;
    double s = lo + random.rand() * (hi - lo);

    int idInteraction;

    if (ID == IDbkg) {
        if (ID > 0) {
            idInteraction = 1; // in Rhorry's "nu_e nu_e -> nu_e nu_e"
        } else {
            idInteraction = 3; // in Rhorry's "nu_e~ nu_e~ -> nu_e~ nu_e~"
        }
    } else {
        if (ID > 0) {
            idInteraction = 5; // in Rhorry's "nu_e nu_mu -> nu_e nu_mu"
        } else {
            idInteraction = 8; // in Rhorry's "nu_e~ nu_mu~ -> nu_e~ nu_mu~"
        }
    }

    // sample the cosine of theta13_com
    double costh13_com = this->secondariesDistribution.sample(idInteraction, s, mass, this->neutrinoFieldMass);
    //std::cout << "costh13_com: " << costh13_com << std::endl;

    std::vector<double> energies;
    double z = candidate->getRedshift();
    double targetMinEnergy = (1 + z) * this->neutrinoField->getMinimumNeutrinoEnergy(z);
    double targetMaxEnergy = (1 + z) * this->neutrinoField->getMaximumNeutrinoEnergy(z);
    try {
        // see if this function wants the neutrino masses in J or kg!
        // Use a local instance to avoid shared mutable state across threads.
        RelativisticInteraction relInteraction(mass / c_squared,
                                               this->neutrinoFieldMass / c_squared,
                                               E, s);
        // energies of the secondary particles
        energies = relInteraction.getProductEnergiesLab(s, costh13_com, mass / c_squared, this->neutrinoFieldMass / c_squared);
        if (energies.size() < 2) {
            std::cout << "No secondaries: energies.size()=" << energies.size()
                      << " channel=" << idInteraction << std::endl;
            return;
        }

        double maxAvailable = E + relInteraction.getTargetEnergyLab();
        double totalOutgoing = energies[0] + energies[1];
        bool invalidEnergy = (!std::isfinite(totalOutgoing) || !std::isfinite(maxAvailable));
        bool highEnergy = (totalOutgoing > 1.05 * maxAvailable);
        for (double e : energies) {
            if (!std::isfinite(e) || e < 0 || e > 1.05 * maxAvailable)
                highEnergy = true;
            if (!std::isfinite(e) || e < 0)
                invalidEnergy = true;
        }
        if (invalidEnergy || highEnergy) {
            std::cout << "Bad secondary energy: E=" << E / EeV
                      << " target=" << relInteraction.getTargetEnergyLab() / EeV
                      << " out=" << totalOutgoing / EeV
                      << " ratio=" << totalOutgoing / maxAvailable
                      << std::endl;
            if (invalidEnergy)
                throw std::runtime_error("Invalid secondary energy");
        }
    } catch (const std::exception& ex) {
        std::cout << "NuNu secondary kinematics failed: " << ex.what()
                  << " E=" << E / EeV
                  << " s=" << s / (GeV * GeV)
                  << " channel=" << idInteraction
                  << " ID=" << ID
                  << " IDbkg=" << IDbkg
                  << " mass=" << mass / eV
                  << " fieldMass=" << this->neutrinoFieldMass / eV
                  << " targetMin=" << targetMinEnergy / eV
                  << " targetMax=" << targetMaxEnergy / eV
                  << " costh13_com=" << costh13_com
                  << std::endl;
        throw;
    }
    if (energies.size() < 2) {
        std::cout << "No secondaries: energies.size()=" << energies.size()
                  << " channel=" << idInteraction << std::endl;
        return;
    }
    for (double e : energies) {
        if (!std::isfinite(e))
            return;
    }

    //std::cout << "Energy0 (EeV)" << energies[0] / EeV << std::endl;
    //std::cout << "Energy1 (EeV)" << energies[1] / EeV << std::endl;

    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

    candidate->current.setEnergy(energies[0]);

    // add the thinning
    if (haveSecondaries) {
        candidate->addSecondary(IDbkg, energies[1], pos, w, interactionTag);
    }
}

double NeutrinoNeutrinoInteraction::findClosestRedshift(double z, const std::vector<double> &redshifts) const {
    auto it = std::lower_bound(redshifts.begin(), redshifts.end(), z);

    if (it == redshifts.begin()) return *it;
    if (it == redshifts.end()) return redshifts.back();

    double upper = *it;
    double lower = *(it - 1);

    return (std::abs(upper - z) < std::abs(lower - z)) ? upper : lower;
}

int NeutrinoNeutrinoInteraction::interactionIndex(int ID, int IDbkg, double mass, double z) const {

    int massIndex = this->neutrinoMixing->massToIndexMass(mass / eV);
    int indexMass = massIndex + 1;
    std::string massComb = this->neutrinoField->getFieldName() + "_m" + std::to_string(indexMass);

    std::string alphaBeta;
    if (ID == IDbkg) {
        alphaBeta = "_NuNu";
    } else {
        alphaBeta = "_NuiNuj";
    }

    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50}; // same as computed in NuPropa-data

    double zClosest = findClosestRedshift(z, redshifts);
    std::ostringstream out;
    out << std::fixed << std::setprecision(1) << zClosest;
    std::string zDec = out.str();
    std::string redshift = "_z" + zDec;

    std::string key = massComb + redshift + alphaBeta;
    auto it = this->ratesDictionary.find(key);

    int index;
    if (it != this->ratesDictionary.end()) {
        index = it->second;
    } else {
        std::ostringstream err;
        err << "NeutrinoNeutrinoInteraction rate key not found: " << key
            << " ID=" << ID
            << " IDbkg=" << IDbkg
            << " mass=" << mass / eV
            << " massIndex=" << massIndex
            << " field=" << this->neutrinoField->getFieldName()
            << " z=" << z
            << " zClosest=" << zClosest;
        throw std::runtime_error(err.str());
    }

    return index;
}

void NeutrinoNeutrinoInteraction::process(Candidate *candidate) const {

    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy();
    int ID = candidate->current.getId();

    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;

    Random &random = Random::instance();
    int sign = (random.rand() < 0.5) ? -1 : +1;
    int IDbkg = sign * this->neutrinoMixing->fromMassToFlavour(this->neutrinoFieldMass / eV);

    //std::cout << "IDbkg: " << IDbkg  << ", mass (eV): " << this->neutrinoFieldMass / eV << std::endl;

    if (!(abs(IDbkg) == 12 || abs(IDbkg) == 14 || abs(IDbkg) == 16))
        return;

    if (!(ID * IDbkg > 0))
        return;

    double mass = this->neutrinoMixing->fromFlavourToMass(ID) * eV; // returned in eV from the function

    std::vector<double> vecEnergy;
    std::vector<double> vecRate;

    int index = interactionIndex(ID, IDbkg, mass, z);

    vecEnergy = this->tabEnergy[index];
    vecRate = this->tabRate[index];

    // check if in tabulated energy range
    if (E < vecEnergy.front() or (E > vecEnergy.back()))
        return;

    // interaction rate
    double rate = interpolate(E, vecEnergy, vecRate);
    if (!std::isfinite(rate) || rate <= 0)
        return;

    // check for interaction
    double randDistance = -log(random.rand()) / rate;
    double step = candidate->getCurrentStep();
    if (step < randDistance) {
        candidate->limitNextStep(limit / rate);
        return;
    } else { // after performing interaction neutrino ceases to exist (hence return)
        performInteraction(candidate, index, mass, IDbkg);
        return;
    }

}

void NeutrinoNeutrinoInteraction::setInteractionTag(std::string tag) {
    this->interactionTag = tag;
}

std::string NeutrinoNeutrinoInteraction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace nupropa
