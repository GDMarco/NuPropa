#ifndef NUPROPA_NEUTRINOPHOTONINTERACTION_H
#define NUPROPA_NEUTRINOPHOTONINTERACTION_H

#include <crpropa/Units.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Module.h>
#include <crpropa/Candidate.h>
#include <crpropa/PhotonBackground.h>

#include <fstream>
#include <cmath>
#include <unordered_map>

namespace nupropa {
/// A custom C++ module for Neutrino-Photon Interaction in astrophysical scenarios. On-shell production of the W boson.
///
using namespace crpropa;

class NeutrinoPhotonInteraction : public Module {
private:
    ref_ptr<PhotonField> photonField;
    bool haveSecondaries;
    double limit;
    //double thinning;
    
    std::string interactionTag = "NGI";
    
    std::vector<std::vector<double>> tabEnergy; //!< neutrino energy in [J]
    std::vector<std::vector<double>> tabRate; //!< interaction rate in [1/m]
    
    std::vector<std::vector<double>> tabE; //!< neutrino energy in [J]
    std::vector<std::vector<double>> tabs; //!< s_kin = s - m^2 in [J**2**]
    std::vector<std::vector<std::vector<double>>> tabCDF; //!< cumulative interaction rate
    
    std::unordered_map<int, std::string> dictionaryNeutrino;
    
public:
    /// The parent's constructor need to be called on initialization!
    NeutrinoPhotonInteraction(ref_ptr<PhotonField> photonField, bool haveSecondaries = false, double limit = 0.1); //double thinning = 0,
    
    // set the target photon field
    void setPhotonField(ref_ptr<PhotonField> photonField);

    // decide if secondary electrons are added to the simulation
    void setHaveSecondaries(bool haveSecondaries);

    /** Limit the propagation step to a fraction of the mean free path
     * @param limit fraction of the mean free path
     */
    void setLimit(double limit);

    
    /** Apply thinning with a given thinning factor
     * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
     */
    /**
    void setThinning(double thinning);
    */
    
    /** set a custom interaction tag to trace back this interaction
     * @param tag string that will be added to the candidate and output
     */
    void setInteractionTag(std::string tag);
    std::string getInteractionTag() const;
    
    void initRate(std::string fname);
    
    std::vector<double> getTabulatedEnergy(int ID) const;
    std::vector<double> getTabulatedRate(int ID) const;
    
    void process(Candidate *candidate) const;
    void performInteraction(Candidate *candidate) const;
};

} // end namespace nupropa

#endif
