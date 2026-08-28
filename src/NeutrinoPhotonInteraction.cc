#include "nupropa/NeutrinoPhotonInteraction.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/ParticleData.h"
#include "nupropa/NeutrinoMixing.h"

#include <crpropa/Units.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Module.h>
#include <crpropa/Candidate.h>
#include <crpropa/PhotonBackground.h>

#include <string>
#include <iomanip>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <filesystem>
#include <unordered_map>
#include <cmath>
#include <iostream>

namespace nupropa {

using namespace crpropa;

NeutrinoPhotonInteraction::NeutrinoPhotonInteraction(ref_ptr<PhotonField> photonField, ref_ptr<NeutrinoMixing> neutrinoMixing, bool haveSecondaries, double limit) : Module(), secondariesDistribution("dsigdcosth") { //double thinning,

    setHaveSecondaries(haveSecondaries);
    setPhotonField(photonField);
    setLimit(limit);
    //setThinning(thinning);
    setNeutrinoMixing(neutrinoMixing);

}

void NeutrinoPhotonInteraction::setPhotonField(ref_ptr<PhotonField> photonField) {

    this->photonField = photonField;
    std::string fname = photonField->getFieldName();
    setDescription("NeutrinoPhotonInteraction::Module" + fname);

    std::string pathPh = "/sdf/home/g/gaetano/CRPropa/CRPropa3-data/dataOff/NeutrinoInteractions/NeutrinoPhotonInteraction/";

    initRate(pathPh);
    initCumulativeRate(pathPh);
    // Disabled for the on-shell W approximation in performInteraction().
    // buildSecondariesDistributionClass();

}

// set the neutrino mixing parameters
void NeutrinoPhotonInteraction::setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    this->neutrinoMixing = neutrinoMixing;
}

void NeutrinoPhotonInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoPhotonInteraction::setLimit(double limit) {
    this->limit = limit;
}

/**
void NeutrinoPhotonInteraction::setThinning(double thinning) {
    this->thinning = thinning;
}
*/

void NeutrinoPhotonInteraction::loadRateFile(const std::string& fileName) {
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

void NeutrinoPhotonInteraction::initRate(std::string filePath) {

    tabEnergy.clear();
    tabRate.clear();
    ratesDictionary.clear();

    std::vector<std::string> flavours = {"Electron", "Muon", "Tau"};
    std::vector<std::string> masses = {"m1", "m2", "m3"};

    int i = 0;
    std::unordered_map<std::string, int> ratesDict;

    // 3 masses for each redshift and each flavour! m1, m2, m1, ...
    /**
     Electron .. _m1_z0     0
     Muon ..  m1_z0     1
     Tauon .. m2_z0     2
    ...
    pathPh = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoPhotonInteraction/"; "NeutrinoElectronPhotonInteraction"  " rate_CMB_m1_z0.0.txt"
    */
    for (const auto& f : flavours) {
        for (const auto& a : masses) {

            loadRateFile(filePath + "Neutrino" + f + "PhotonInteraction/rate_" + this->photonField->getFieldName() + "_" +  a + ".txt");
            ratesDict[f + "_" + this->photonField->getFieldName() + "_" + a] = i;

            i = i + 1;

        }
    }

    this->ratesDictionary = ratesDict;

}

void NeutrinoPhotonInteraction::loadCumulativeRateFile(const std::string& fileName) {

    std::ifstream infile(fileName.c_str());

    if (!infile.good())
        throw std::runtime_error("NeutrinoPhotonInteraction: could not open file" + fileName);

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

void NeutrinoPhotonInteraction::initCumulativeRate(std::string filePath) {

    tabE.clear();
    tabs.clear();
    tabCDF.clear();

    std::vector<std::string> flavours = {"Electron", "Muon", "Tau"};
    std::vector<std::string> masses = {"m1", "m2", "m3"};

    // 3 masses for each redshift and each flavour! m1, m2, m1, ...
    /**
     Electron .. _m1_z0     0
     Muon ..  m1_z0     1
     Tauon .. m2_z0     2
    ...
    pathPh = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoPhotonInteraction/"; "NeutrinoElectronPhotonInteraction"  " cdf_CMB_m1_z0.0.txt"
    */
    for (const auto& f : flavours) {
        for (const auto& a : masses) {

            loadCumulativeRateFile(filePath + "Neutrino" + f + "PhotonInteraction/cdf_" + this->photonField->getFieldName() + "_" +  a + ".txt");

        }
    }
}

void NeutrinoPhotonInteraction::buildSecondariesDistributionClass() {

    if (not haveSecondaries)
        return;

    /**if (ID == 12 || ID == 14 || ID == 16)
     return (ID + 44) / 2;   // neutrinos: 12 to 28, 14 to 29, 16 to 30
 if (ID == -12 || ID == -14 || ID == -16)
     return (50 - ID) / 2;   // antineutrinos: -12 to 31, -14 to 32, -16 to 33
 return -1; // invalid ID*/

    std::vector<int> channelsId = {28, 29, 30, 31, 32, 33};

    for (int i = 0; i < channelsId.size(); ++i) {
        this->secondariesDistribution.buildChannel(channelsId[i]);
    }
}

void NeutrinoPhotonInteraction::setRelativisticInteraction(double m1, double E, double s) const {
    this->relInteraction = new RelativisticInteraction(m1, E, s);
}

int NeutrinoPhotonInteraction::interactionIndex(int ID, double mass) const {

    int indexMass = this->neutrinoMixing->massToIndexMass(mass / eV) + 1;
    std::string massComb = this->photonField->getFieldName() + "_m" + std::to_string(indexMass);

    std::string alpha;
    if (abs(ID) == 12) {
        alpha = "Electron";
    } else if (abs(ID) == 14) {
        alpha = "Muon";
    } else if (abs(ID) == 16) {
        alpha = "Tau";
    } else {
        throw std::runtime_error("Unsupported neutrino ID in interactionIndex: " + std::to_string(ID));
    }

    auto it = this->ratesDictionary.find(alpha + "_" + massComb);

    int index;
    if (it != this->ratesDictionary.end()) {
        index = it->second;
    } else {
        throw std::runtime_error("Index not found in the dictionary rate!");
    }

    return index;
}

int NeutrinoPhotonInteraction::fromIDtoChannel(int ID) const {
    if (ID == 12 || ID == 14 || ID == 16)
        return (ID + 44) / 2;   // neutrinos: 12 to 28, 14 to 29, 16 to 30
    if (ID == -12 || ID == -14 || ID == -16)
        return (50 - ID) / 2;   // antineutrinos: -12 to 31, -14 to 32, -16 to 33
    return -1; // invalid ID
}

void NeutrinoPhotonInteraction::performInteraction(Candidate *candidate, int index, double mass) const {

    if (not haveSecondaries) {
        candidate->setActive(false);
        return;
    }

    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    int ID = candidate->current.getId();
    double w = 1.; // no weights so far!

    Random &random = Random::instance();
    int WbosonID = (ID > 0) ? 24 : -24;
    int leptonID = (ID > 0) ? ID - 1 : ID + 1;
    ParticleData particle;
    double WbosonRestEnergy = particle.getParticleMass(std::abs(WbosonID)) * c_squared;
    double leptonRestEnergy = particle.getParticleMass(std::abs(leptonID)) * c_squared;
    // Work in the interaction frame at redshift z. The W is produced exactly
    // on shell and handed to the decay plugin without redshift rescaling; the
    // charged lepton receives the remaining interaction energy and follows the
    // usual CRPropa storage convention.
    const double pythiaThresholdMargin = 1e-6;
    double WbosonEnergy = WbosonRestEnergy;
    double leptonEnergy = E - WbosonEnergy;
    Vector3d pos = random.randomInterpolatedPosition(
        candidate->previous.getPosition(), candidate->current.getPosition());

    candidate->setActive(false);
    // Produce the W on shell and give the exact remaining interaction energy
    // to the corresponding charged lepton.
    candidate->addSecondary(WbosonID, WbosonEnergy, pos, w, interactionTag);
    if (leptonEnergy < leptonRestEnergy * (1 + pythiaThresholdMargin))
        return;
    candidate->addSecondary(leptonID, leptonEnergy / (1 + z), pos, w, interactionTag);

#if 0
    // Differential-cross-section secondary production, temporarily disabled.
    std::vector<double> vecE = this->tabE[index];
    std::vector<double> vecs = this->tabs[index];
    std::vector<std::vector<double>> vecCDF = this->tabCDF[index];

    // check if in tabulated energy range
    if (vecE.empty() || vecs.empty() || vecCDF.empty()) {
        std::cout << "NuPhoton no secondaries: empty tables"
                  << " index=" << index << std::endl;
        return;
    }
    if (E < vecE.front() or (E > vecE.back())) {
        std::cout << "NuPhoton no secondaries: E outside secondary table"
                  << " E=" << E / EeV
                  << " Emin=" << vecE.front() / EeV
                  << " Emax=" << vecE.back() / EeV
                  << std::endl;
        return;
    }

    int leptonID;
    int WbosonID;

    int channelID = fromIDtoChannel(ID);

    if (ID > 0) {
        leptonID = ID - 1;
        WbosonID = 24; // W+
    } else {
        leptonID = ID + 1;
        WbosonID = -24; // W-
    }

    ParticleData particle;
    double leptonMass = particle.getParticleMass(leptonID); // in kg
    double WbosonMass = particle.getParticleMass(WbosonID); // in kg

    double sThr = std::pow(leptonMass + WbosonMass, 2) * c_squared * c_squared;

    // sample the value of s
    size_t i = closestIndex(E, vecE);  // find closest tabulation point
    if (i >= vecCDF.size()) {
        std::cout << "NuPhoton no secondaries: CDF index out of range"
                  << " i=" << i << " size=" << vecCDF.size() << std::endl;
        return;
    }
    size_t j = random.randBin(vecCDF[i]);
    if (j >= vecs.size()) {
        std::cout << "NuPhoton no secondaries: s-bin out of range"
                  << " j=" << j << " size=" << vecs.size() << std::endl;
        return;
    }

    double hi = vecs[j];
    double lo = (j == 0) ? sThr : std::max(sThr, vecs[j - 1]); // first s-tabulation point below min(s_kin)
    if (hi < lo)
        hi = lo;
    double s = lo + random.rand() * (hi - lo); // should I add the neutrino masses? since it is the s_kin!!

    // sample the cosine of theta13_com
    double costh13_com = this->secondariesDistribution.sample(channelID, s, sThr);

    std::vector<double> energies;
    double targetMinEnergy = photonField->getMinimumPhotonEnergy(z);
    double targetMaxEnergy = photonField->getMaximumPhotonEnergy(z);
    try {
        // Use a local instance to avoid shared mutable state across threads.
        RelativisticInteraction relInteraction(mass / c_squared, E, s);
        // energies of the secondary particles
        energies = relInteraction.getProductEnergiesLab(s, costh13_com, leptonMass, WbosonMass);
        if (energies.size() < 2) {
            std::cout << "No secondaries: energies.size()=" << energies.size()
                      << " channel=" << channelID << std::endl;
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
        std::cout << "NuPhoton secondary kinematics failed: " << ex.what()
                  << " E=" << E / EeV
                  << " s=" << s / (GeV * GeV)
                  << " channel=" << channelID
                  << " leptonID=" << leptonID
                  << " WbosonID=" << WbosonID
                  << " leptonMass=" << leptonMass
                  << " WbosonMass=" << WbosonMass
                  << " targetMin=" << targetMinEnergy / eV
                  << " targetMax=" << targetMaxEnergy / eV
                  << " costh13_com=" << costh13_com
                  << std::endl;
        throw;
    }
    if (energies.size() < 2) {
        std::cout << "No secondaries: energies.size()=" << energies.size()
                  << " channel=" << channelID << std::endl;
        return;
    }
    for (double e : energies) {
        if (!std::isfinite(e))
            return;
    }

    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->setActive(false);

    //if (haveSecondaries)?, add thinning
    candidate->addSecondary(leptonID, energies[0] / (1 + z), pos, w, interactionTag);
    // W production, to check if on shell?
    candidate->addSecondary(WbosonID, energies[1] / (1 + z), pos, w, interactionTag);
#endif

}

void NeutrinoPhotonInteraction::process(Candidate *candidate) const {

    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    double ID = candidate->current.getId();

    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;

     double mass = this->neutrinoMixing->fromFlavourToMass(ID) * eV; // returned in eV from the function

    int index = interactionIndex(ID, mass);

    std::vector<double> vecEnergy = this->tabEnergy[index];
    std::vector<double> vecRate = this->tabRate[index];

    // check if in tabulated energy range
    if (E < vecEnergy.front() or (E > vecEnergy.back()))
        return;

    // interaction rate
    double rate = interpolate(E, vecEnergy, vecRate);
    if (!std::isfinite(rate) || rate <= 0)
        return;
    rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

    // check for interaction
    Random &random = Random::instance();
    double randDistance = -log(random.rand()) / rate;
    double step = candidate->getCurrentStep();
    if (step < randDistance) {
        candidate->limitNextStep(limit / rate);
        return;
    } else { // after performing interaction neutrino ceases to exist (hence return)
        performInteraction(candidate, index, mass);
        return;
    }
}

void NeutrinoPhotonInteraction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string NeutrinoPhotonInteraction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace nupropa
