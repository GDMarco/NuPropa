#ifndef NUPROPA_HADRONISATIONS_H
#define NUPROPA_HADRONISATIONS_H

#include <crpropa/Module.h>
#include <crpropa/Candidate.h>

#include<string>

namespace nupropa {

using namespace crpropa;

class Hadronisations : public Module {
private:
    
    bool haveOtherSecondaries;
    bool haveNeutrinos;
    bool angularCorrection;
    double limit;
    // double thinning;
    mutable std::string hadronisationTag = "Hadr";
    mutable bool hasPending = false;
    mutable std::vector<int> pendingIds;
    mutable std::vector<double> pendingEs;
    mutable Vector3d pendingPos;
    mutable double pendingW = 1.0;

public:
    
    /** Constructor
     @param haveOtherSecondaries    if true, add secondary particles as candidates that are not neutrinos
     @param haveNeutrinos    if true, add secondary neutrinos as candidates
     @param angularCorrection    if true, consider the decay angular distribution of the secondaries
     // @param thinning        weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
     @param limit            step size limit as fraction of mean free path
     */
    Hadronisations(bool haveOtherSecondaries = false, bool haveNeutrinos = false, bool angularCorrection = false); // double thinning = 0,
    
    void setHaveOtherSecondaries(bool haveOtherSecondaries);
    
    void setHaveNeutrinos(bool haveNeutrinos);
    
    void setAngularCorrection(bool angularCorrection);
    
    /** Apply thinning with a given thinning factor
     * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
     */
    // void setThinning(double thinning);
    
    void setHadronisationTag(std::string tag) const;
    std::string getHadronisationTag() const;
    
    void performHadronisation(Candidate *candidate, std::vector<int> Ids, std::vector<double> Es, Vector3d pos, double w) const;
    void setPendingHadronisation(std::vector<int> Ids, std::vector<double> Es, Vector3d pos, double w) const;
    void process(Candidate *candidate) const override;
    
};

} // end namespace nupropa
#endif
