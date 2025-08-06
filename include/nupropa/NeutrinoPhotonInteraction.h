#ifndef NUPROPA_NEUTRINOPHOTONINTERACTION_H
#define NUPROPA_NEUTRINOPHOTONINTERACTION_H

#include <crpropa/Units.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Module.h>
#include <crpropa/Candidate.h>
#include <crpropa/PhotonBackground.h>

#include "nupropa/NeutrinoMixing.h"
#include "nupropa/RelativisticInteraction.h"

#include <fstream>
#include <cmath>
#include <unordered_map>

namespace nupropa {

using namespace crpropa;

class NeutrinoPhotonInteraction : public Module {
private:
    
    ref_ptr<PhotonField> photonField;
    ref_ptr<NeutrinoMixing> neutrinoMixing;
    
    bool haveSecondaries;
    double limit;
    // double thinning;
    
    std::string interactionTag = "NuPhI";
    
    std::vector<std::vector<double>> tabEnergy; //!< neutrino energy in [J]
    std::vector<std::vector<double>> tabRate; //!< interaction rate in [1/m]
    std::vector<std::vector<double>> tabE; //!< neutrino energy in [J]
    std::vector<std::vector<double>> tabs; //!< s_kin = s - m^2 in [J**2**]
    std::vector<std::vector<std::vector<double>>> tabCDF; //!< cumulative interaction rate
    std::unordered_map<std::string, int> ratesDictionary;
    
    mutable ref_ptr<RelativisticInteraction> relInteraction;
    
public:
    
    NeutrinoPhotonInteraction(ref_ptr<PhotonField> photonField, ref_ptr<NeutrinoMixing> neutrinoMixing, bool haveSecondaries = false, double limit = 0.1); //double thinning = 0,
    
    // set the target photon field
    void setPhotonField(ref_ptr<PhotonField> photonField);

    // set the neutrino mixing parameters
    void setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing);
    
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
    
    void loadRateFile(const std::string& fileName);
    void loadCumulativeRateFile(const std::string& fileName);
    
    void initRate(std::string filePath);
    void initCumulativeRate(std::string filePath);
    
    void setRelativisticInteraction(double m1, double E, double s) const;
    
    int fromIDtoChannel(int ID) const;
    int interactionIndex(int ID, double mass) const;
    
    void process(Candidate *candidate) const;
    void performInteraction(Candidate *candidate, int index, double mass) const;
};

} // end namespace nupropa

#endif
