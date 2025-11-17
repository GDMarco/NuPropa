#include "nupropa/NeutrinoOscillation.h"
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

NeutrinoOscillation::NeutrinoOscillation() {};

NeutrinoOscillation::NeutrinoOscillation(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    setNeutrinoMixing(neutrinoMixing);
}

void NeutrinoOscillation::setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    this->neutrinoMixing = neutrinoMixing;
}

void NeutrinoOscillation::process(Candidate *candidate) const {
    
    int ID = candidate->current.getId();
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;
    
    double E = candidate->current.getEnergy();
    double step = candidate->getCurrentStep();
    // it should be taken as the average probabilities!

    Vector3d pos1 = candidate->previous.getPosition();
    Vector3d pos2 = candidate->current.getPosition();
    
    int newID = this->neutrinoMixing->oscillateFlavour(ID, E, step);

    // candidate->limitNextStep(1 * kpc);
    if (newID != ID) {
        // Optionally, randomize the position along the step
        Random &random = Random::instance();
        Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
        candidate->current.setPosition(pos);
        candidate->current.setId(newID);
        std::cout << "Oscillated to ID: " << newID << std::endl;
    }
    //evolve(candidate, ID, E, step);
    // Formal note:
    // current is a ParticleState object. To take it properly I need to adjust how they associate the mass to the ID for neutrinos! (Differently wrt the other particles.)
    // For neutrinos, the mixingNeutrino class should be used to determine the masses, directly in the particleState!
}


} // end namespace nupropa

