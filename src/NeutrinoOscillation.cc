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

void NeutrinoOscillation::process(crpropa::Candidate *candidate) const {
    
    int ID = candidate->current.getId();
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;
    
    // double z = candidate->getRedshift();
    double E = candidate->current.getEnergy();
    double L = candidate->getTrajectoryLength();
    
    int IDosc = this->neutrinoMixing->oscillateFlavour(ID, E, L);
    
    candidate->current.setId(IDosc);
    // Formal note: 
    // current is a ParticleState object. To take it properly I need to adjust how they associate the mass to the ID for neutrinos! (Differently wrt the other particles.)
    // For neutrinos, the mixingNeutrino class should be used to determine the masses, directly in the particleState!
    
}


} // end namespace nupropa

