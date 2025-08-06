#ifndef NUPROPA_NEUTRINOOSCILLATION_H
#define NUPROPA_NEUTRINOOSCILLATION_H

#include "nupropa/NeutrinoMixing.h"

#include <crpropa/Module.h>
#include <crpropa/Candidate.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>


#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace nupropa {

using namespace crpropa;

class NeutrinoOscillation : public Module
{
private:
    
    ref_ptr<NeutrinoMixing> neutrinoMixing;
    
public:
    
    NeutrinoOscillation();
    NeutrinoOscillation(ref_ptr<NeutrinoMixing> neutrinoMixing);
    
    void setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing);
    
    void process(crpropa::Candidate *candidate) const;
    
};

} // end namespace nupropa

#endif
