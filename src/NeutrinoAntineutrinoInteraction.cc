#include "nupropa/NeutrinoAntineutrinoInteraction.h"
#include "nupropa/Channels.h"
#include "nupropa/ChannelsBundle.h"
#include "nupropa/NeutrinoField.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/ParticleData.h"
#include "nupropa/NeutrinoMixing.h"
#include "nupropa/Hadronisations.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include <iomanip>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace nupropa {

using namespace crpropa;

NeutrinoAntineutrinoInteraction::NeutrinoAntineutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, ref_ptr<Channels> channels, ref_ptr<NeutrinoMixing> neutrinoMixing, bool haveSecondaries, double limit) : Module(), secondariesDistribution("dsigdcosth") { // double thinning

    setChannels(channels);
    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    //setThinning(thinning);
    setNeutrinoField(neutrinoField);
    setNeutrinoMixing(neutrinoMixing);

}

/**
 One assumes that the tables have been produced within the same computation, so the energy bins are the same.
 In the case of the neutrino-antineutrino interactions the elastic scatterings furnish the tabEnergy, since they do not have interaction energy thresholds (almost)!
 */

void NeutrinoAntineutrinoInteraction::setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField) {
    this->neutrinoField = neutrinoField;
    this->neutrinoFieldMass = neutrinoField->getMass();
    std::string fname = neutrinoField->getFieldName();
    setDescription("NeutrinoAntineutrinoInteraction::Module" + fname);
    setInteractionTag("NuAntiNuInt");

    setChannelsBundle(this->channels, fname);

    buildSecondariesDistributionClass();
}

void NeutrinoAntineutrinoInteraction::setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    this->neutrinoMixing = neutrinoMixing;
}

void NeutrinoAntineutrinoInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoAntineutrinoInteraction::setLimit(double limit) {
    this->limit = limit;
}

/**
void NeutrinoNeutrinoInteraction::setThinning(double thinning) {
    this->thinning = thinning;
}
*/

 void NeutrinoAntineutrinoInteraction::setChannels(ref_ptr<Channels> channels) {
     this->channels = channels;
}

void NeutrinoAntineutrinoInteraction::setChannelsBundle(ref_ptr<Channels> channels, std::string fname) {
    this->channelsBundle = new ChannelsBundle(channels, fname);
}

void NeutrinoAntineutrinoInteraction::setRelativisticInteraction(double m1, double m2, double E, double s) const {
    this->relInteraction = new RelativisticInteraction(m1, m2, E, s);
}

void NeutrinoAntineutrinoInteraction::buildSecondariesDistributionClass() {

    if (!haveSecondaries)
        return;

    std::vector<std::string> interactionChannels = this->channels->getInteractionChannels();
    std::vector<int> activeChannels = this->channels->getActiveChannels();
    std::string path = channels->getInteractionFolderPath();

    for (int i = 0; i < interactionChannels.size(); ++i) {

        if (activeChannels[i] == 0)
            continue;

        std::string pathInt = path + interactionChannels[i];
        std::string filenameProdChan = pathInt + "/products_channelId.txt";

        std::vector<int> idProdChan = getProductsChannelId(filenameProdChan);

        this->secondariesDistribution.buildChannel(idProdChan[2]);

    }
}

std::vector<int> NeutrinoAntineutrinoInteraction::getProductsChannelId(std::string& filename) {

    std::ifstream infile(filename.c_str());

    if (!infile.good())
        throw std::runtime_error("Could not open rate file: " + filename);

    std::vector<int> Ids;

    while (infile.good()) {
        if (infile.peek() != '#') {
            int a, b, c;
            infile >> a >> b >> c; // 2 to 2 process, only one row (to optimise)
            if (infile) {
                Ids = {a, b, c};
            }
        }
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infile.close();

    return Ids;
}

bool NeutrinoAntineutrinoInteraction::isQuark(int pdgId) const {
    int id = std::abs(pdgId);
    return (id >= 1 && id <= 6);
}

void NeutrinoAntineutrinoInteraction::performInteraction(Candidate *candidate, double mass, int IDbkg,
                                                         const ChannelsBundle::ChannelSelection &selection) const {

    double E = candidate->current.getEnergy();
    int ID = candidate->current.getId();

    if (not haveSecondaries) {
        candidate->setActive(false);
        return;
    }

    //int indexChannel = this->channelsBundle->getSelectedIndex();

    //std::cout << "indexChannel: " << indexChannel << std::endl;

    const std::vector<double> &tabE = selection.tabE;
    const std::vector<double> &tabs = selection.tabs;
    const std::vector<std::vector<double>> &tabCDF = selection.tabCDF;
    const std::vector<int> &prodChanId = selection.prodChanId; // the channel ID is inside the products.txt file in each folder.

    // check if in tabulated energy range
    if (tabE.empty() || tabs.empty() || tabCDF.empty()) {
        std::cout << "NuAntiNu no secondaries: empty tables"
                  << " channel=" << (prodChanId.size() > 2 ? prodChanId[2] : -1)
                  << std::endl;
        return;
    }
    if (E < tabE.front() or (E > tabE.back())) {
        std::cout << "NuAntiNu no secondaries: E outside secondary table"
                  << " E=" << E / EeV
                  << " Emin=" << tabE.front() / EeV
                  << " Emax=" << tabE.back() / EeV
                  << std::endl;
        return;
    }

    int prodId1 = prodChanId[0];
    int prodId2 = prodChanId[1];
    int chanId = prodChanId[2];

    // for the neutrino antineutrino Z resonance into nu nux of another flavor
    Random &random = Random::instance();
    if (chanId == 9) {
        if (std::abs(ID) == 12) {
            int Id = (random.rand() < 0.5) ? 14 : 16;
            prodId1 = Id;
            prodId2 = -Id;
        } else if (std::abs(ID) == 14) {
            int Id = (random.rand() < 0.5) ? 12 : 16;
            prodId1 = Id;
            prodId2 = -Id;
        } else { // abs(ID) = 16
            int Id = (random.rand() < 0.5) ? 14 : 12;
            prodId1 = Id;
            prodId2 = -Id;
        }
    }

    double m3;
    double m4;

    if (prodId1 == 0 && prodId2 == 0) {  // in the elastic scattering they are set to 0!

        prodId1 = prodId1 + ID;
        prodId2 = prodId2 + IDbkg;

        m3 = mass / c_squared; // it is given in J
        m4 = this->neutrinoFieldMass / c_squared; // it is given in J

    } else {

        ParticleData particle;
        auto productMass = [&](int pdgId) {
            int absId = std::abs(pdgId);
            if (absId == 12 || absId == 14 || absId == 16)
                return this->neutrinoMixing->fromFlavourToMass(pdgId) * eV / c_squared;
            return particle.getParticleMass(absId); // in kg
        };
        m3 = productMass(prodId1);
        m4 = productMass(prodId2);
    }

    double sThr = std::max(std::pow(m3 + m4, 2) * c_squared * c_squared, 1e-100);

    // sample the value of s
    size_t i = closestIndex(E, tabE);  // find closest tabulation point
    if (i >= tabCDF.size()) {
        std::cout << "NuAntiNu no secondaries: CDF index out of range"
                  << " i=" << i << " size=" << tabCDF.size() << std::endl;
        return;
    }
    size_t j = random.randBin(tabCDF[i]);
    if (j >= tabs.size()) {
        std::cout << "NuAntiNu no secondaries: s-bin out of range"
                  << " j=" << j << " size=" << tabs.size() << std::endl;
        return;
    }
    double hi = tabs[j];
    double lo = (j == 0) ? sThr : std::max(sThr, tabs[j - 1]); // first s-tabulation point below min(s_kin)
    if (hi < lo)
        hi = lo;
    double s = lo + random.rand() * (hi - lo);

    // sample the cosine of theta13_com
    std::string variable = "dsigdcosth";

    // sample the cosine of theta13_com
    double costh13_com = this->secondariesDistribution.sample(chanId, s, sThr);

    std::vector<double> energies;
    double z = candidate->getRedshift();
    double targetMinEnergy = (1 + z) * this->neutrinoField->getMinimumNeutrinoEnergy(z);
    double targetMaxEnergy = (1 + z) * this->neutrinoField->getMaximumNeutrinoEnergy(z);
    try {
        // Use a local instance to avoid shared mutable state across threads.
        RelativisticInteraction relInteraction(mass / c_squared,
                                               this->neutrinoFieldMass / c_squared,
                                               E, s);
        // energies of the secondary particles
        energies = relInteraction.getProductEnergiesLab(s, costh13_com, m3, m4);
        if (energies.size() < 2) {
            std::cout << "No secondaries: energies.size()=" << energies.size()
                      << " channel=" << chanId << std::endl;
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
        std::cout << "NuAntiNu secondary kinematics failed: " << ex.what()
                  << " E=" << E / EeV
                  << " s=" << s / (GeV * GeV)
                  << " channel=" << chanId
                  << " prodId1=" << prodId1
                  << " prodId2=" << prodId2
                  << " m3=" << m3
                  << " m4=" << m4
                  << " targetMin=" << targetMinEnergy / eV
                  << " targetMax=" << targetMaxEnergy / eV
                  << " costh13_com=" << costh13_com
                  << std::endl;
        throw;
    }
    if (energies.size() < 2) {
        std::cout << "No secondaries: energies.size()=" << energies.size()
                  << " channel=" << chanId << std::endl;
        return;
    }
    for (double e : energies) {
        if (!std::isfinite(e)) {
            return; // skip invalid kinematics
        }
    }

    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    candidate->setActive(false);

    // thinning not implemented!
    double w = 1.;
    if (isQuark(prodId1) && isQuark(prodId2)) {
        std::vector<int> Ids = {prodId1, prodId2};
        Hadronisations hadronisations;
        hadronisations.setHaveOtherSecondaries(haveSecondaries);
        hadronisations.setHaveNeutrinos(haveSecondaries);
        hadronisations.setPendingHadronisation(Ids, energies, pos, w);
        hadronisations.process(candidate);
    } else if ((isQuark(prodId1) && !isQuark(prodId2)) ||
               (!isQuark(prodId1) && isQuark(prodId2)) ) {
        throw std::runtime_error("No qqbar state! Only one quark produced.");
    } else {
        candidate->addSecondary(prodId1, energies[0], pos, w, this->interactionTag);
        candidate->addSecondary(prodId2, energies[1], pos, w, this->interactionTag);
    }
}

void NeutrinoAntineutrinoInteraction::process(Candidate *candidate) const {

    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy();
    int ID = candidate->current.getId();

    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;

    Random &random = Random::instance();
    int sign = (random.rand() < 0.5) ? -1 : +1;
    int IDbkg = sign * this->neutrinoMixing->fromMassToFlavour(this->neutrinoFieldMass / eV);

    if (ID * IDbkg > 0)
        return;

    double mass = this->neutrinoMixing->fromFlavourToMass(ID) * eV; // returned in eV from the function
    int indexMass = this->neutrinoMixing->massToIndexMass(mass / eV) + 1;

    std::string massComb = this->neutrinoField->getFieldName() + "_m" + std::to_string(indexMass);

    double rate = this->channelsBundle->getRate(ID, IDbkg, massComb, z, E);
    if (!std::isfinite(rate) || rate <= 0)
        return;

    // check for interaction
    double randDistance = -log(random.rand()) / rate;
    double step = candidate->getCurrentStep();
    if (step < randDistance) {
        candidate->limitNextStep(limit / rate);
        return;
    } else { // after performing interaction neutrino ceases to exist (hence return)
        ChannelsBundle::ChannelSelection selection =
            this->channelsBundle->selectInteraction(ID, IDbkg, massComb, z, E);
        if (!selection.valid) {
            std::cout << "NuAntiNu no secondaries: invalid channel selection"
                      << std::endl;
            return;
        }
        performInteraction(candidate, mass, IDbkg, selection);
        return;
    }
}

void NeutrinoAntineutrinoInteraction::setInteractionTag(std::string tag) const {
    this->interactionTag = tag;
}

std::string NeutrinoAntineutrinoInteraction::getInteractionTag() const {
    return this->interactionTag;
}

}
