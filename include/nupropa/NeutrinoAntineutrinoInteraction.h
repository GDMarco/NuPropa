#ifndef NUPROPA_NEUTRINOANTINEUTRINOINTERACTION_H
#define NUPROPA_NEUTRINOANTINEUTRINOINTERACTION_H

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include "nupropa/NeutrinoBackground.h"
#include "nupropa/Channels.h"
#include "nupropa/ChannelsBundle.h"

#include <string>


namespace nupropa {

using namespace crpropa;
/// A custom C++ module for Neutrino-Neutrino Interactions in astorphysical scenarios
class NeutrinoAntineutrinoInteraction : public Module
{
private:
    
    ref_ptr<NeutrinoField> neutrinoField;
    ref_ptr<Channels> channels;
    ref_ptr<ChannelsBundle> channelsBundle;
    bool haveSecondaries;
    double limit;
    //double thinning;
    mutable std::string interactionTag;
    
    int neutrinoFieldID = 12; // check the initialization!
    
    //std::vector<std::vector<double>> tabEnergy; // 4 columns table depending on the neutrino flavour and number, in initRate they should be built properly
    // mutable std::vector<std::vector<double>> tabRate;
    
    // std::vector<std::vector<double>> tabE;
    // std::vector<std::vector<double>> tabs;
    // std::vector<std::vector<std::vector<double>>> tabCDF;
    // mutable std::vector<std::vector<int>> tabProductsID;
    // mutable std::vector<int> selectedProductsID;
    
    // std::unordered_map<int, std::string> interactionDictionary;
    // mutable std::vector<std::vector<double>> channelProbability; // to be sync with interactionDictionary
public:
    /// The parent's constructor need to be called on initialization!
    NeutrinoAntineutrinoInteraction(ref_ptr<NeutrinoField>, ref_ptr<Channels> channels, bool haveSecondaries, double limit); // double thinning = 0,
    
    // set the target neutrino field
    void setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField);
    
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
    
    /** Activate/deactivate the neutrino interaction channels
     * @param channels flags
     */
    void setChannels(ref_ptr<Channels> channels);
    void setChannelsBundle(ref_ptr<Channels> channels, std::string fname);
    
    /** set a custom interaction tag to trace back this interaction
     * @param tag string that will be added to the candidate and output
     */
    void setInteractionTag(std::string tag) const;
    std::string getInteractionTag() const;
    
    void process(crpropa::Candidate *candidate) const;
    void performInteraction(Candidate *candidate) const;
    
};

} // end namespace crpropa
#endif
