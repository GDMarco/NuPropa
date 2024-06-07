#ifndef NEUTRINONEUTRINOINTERACTION_H
#define NEUTRINONEUTRINOINTERACTION_H

#include <crpropa/Module.h>
#include "nupropa/NeutrinoBackground.h"

#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace nupropa {

using namespace crpropa;
/// A custom C++ module for Neutrino-Neutrino Interactions in astorphysical scenarios
class NeutrinoNeutrinoInteraction : public Module
{
private:
    
    ref_ptr<NeutrinoField> neutrinoField;
    int neutrinoFieldID = 12; // check the initialization!
    bool haveSecondaries;
    double limit;
    std::string interactionTag = "NuNuI";
    
    std::vector<std::vector<double>> tabEnergy; // 4 columns table depending on the neutrino flavour and number, in initRate they should be built properly
    std::vector<std::vector<double>> tabRate;
    
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    
public:
    /// The parent's constructor need to be called on initialization!
    NeutrinoNeutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, bool haveSecondaries = false, double limit = 0.1);
    
    // set the target neutrino field
    void setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField);
    
    // decide if secondary neutrinos are added to the simulation
    void setHaveSecondaries(bool haveSecondaries);
    
    /** Limit the propagation step to a fraction of the mean free path
     * @param limit fraction of the mean free path
     */
    void setLimit(double limit);
    
    /** set a custom interaction tag to trace back this interaction
     * @param tag string that will be added to the candidate and output
     */
    void setInteractionTag(std::string tag);
    std::string getInteractionTag() const;
    
    void initRate(std::string fileNuNu, std::string fileNuiNuj);
    // void initCumulativeRate(std::string filename);
    
    void process(crpropa::Candidate *candidate) const;
    void performInteraction(Candidate *candidate) const;
    
};

} // end namespace nupropa

#endif
