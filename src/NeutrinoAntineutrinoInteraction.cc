#include "nupropa/NeutrinoAntineutrinoInteraction.h"
#include "nupropa/Channels.h"
#include "nupropa/ChannelsBundle.h"
#include "nupropa/NeutrinoBackground.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include <string>
#include <filesystem>
#include <unordered_map>
#include <fstream>
#include <vector>

namespace nupropa {

using namespace crpropa;

// The parent's constructor need to be called on initialization!
NeutrinoAntineutrinoInteraction::NeutrinoAntineutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, ref_ptr<Channels> channels, bool haveSecondaries, double limit) : Module() { // double thinning
    setChannels(channels);
    setNeutrinoField(neutrinoField);
    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    //setThinning(thinning);
}

/**
 One assumes that the tables have been produced within the same computation, so the energy bins are the same.
 
 In the case of the neutrino-antineutrino interactions the elastic scatterings furnish the tabEnergy, since they do not have interaction energy thresholds (almost)!
 */

// channels should be an unordered map <bool, "interaction">, maybe interaction tag? setChannels should be a way of filling this object, knowing how they are ordered!
void NeutrinoAntineutrinoInteraction::setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField) {
    this->neutrinoField = neutrinoField;
    this->neutrinoFieldID = neutrinoField->getParticleID(); 
    std::string fname = neutrinoField->getFieldName();
    
    setDescription("NeutrinoAntineutrinoInteraction::Module" + fname);
    
    setChannelsBundle(this->channels, fname);
    // initRate(fname);
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

    this->channelsBundle = new ChannelsBundle(channels, fname); //ChannelsBundle(channels, fname);
}

void NeutrinoAntineutrinoInteraction::performInteraction(Candidate *candidate) const {
    
    candidate->setActive(false);
    
    if (not haveSecondaries)
        return;
    
    // function to take the cumulative inverse mean free path and the other stuff!
    
    // the products depends on the channel
    double w = 1.;
    
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    int ID = candidate->current.getId();
    
    // energies of the secondary particles
    Random &random = Random::instance();
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    // naive values
    double E1 = E / 10.;
    double E2 = E - E1;
    
    if (haveSecondaries)
        candidate->addSecondary(this->channelsBundle->getSelectedProductsID()
                                
                                [0], E1 / (1 + z), pos, w, this->interactionTag);
        candidate->addSecondary(this->channelsBundle->getSelectedProductsID()[1], E2 / (1 + z), pos, w, this->interactionTag);
}

void NeutrinoAntineutrinoInteraction::process(Candidate *candidate) const {
    
    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = (1 + z) * candidate->current.getEnergy();
    int ID = candidate->current.getId();
    int IDBkg = this->neutrinoFieldID;
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16) || !(ID * IDBkg < 0))
        return;

    double rate = this->channelsBundle->getRate(ID, IDBkg, E);
    
    rate *= pow_integer<2>(1 + z) * neutrinoField->getRedshiftScaling(z); // see how it is implemented this function in NeutrinoBackground.h
    
    // check for interaction
    Random &random = Random::instance();
    double randDistance = -log(random.rand()) / rate;
    double step = candidate->getCurrentStep();
    if (step < randDistance) {
        candidate->limitNextStep(limit / rate);
        return;
    } else { // after performing interaction neutrino ceases to exist (hence return)
        performInteraction(candidate);
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

