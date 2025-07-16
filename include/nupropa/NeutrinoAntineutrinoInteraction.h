#ifndef NUPROPA_NEUTRINOANTINEUTRINOINTERACTION_H
#define NUPROPA_NEUTRINOANTINEUTRINOINTERACTION_H

#include <crpropa/Module.h>
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include "nupropa/NeutrinoBackground.h"
#include "nupropa/Channels.h"
#include "nupropa/ChannelsBundle.h"
#include "nupropa/RelativisticInteraction.h"

#include <string>


namespace nupropa {

using namespace crpropa;

class NeutrinoAntineutrinoInteraction : public Module
{
private:
    
    ref_ptr<NeutrinoField> neutrinoField;
    ref_ptr<NeutrinoMixing> neutrinoMixing;
    
    ref_ptr<Channels> channels;
    ref_ptr<ChannelsBundle> channelsBundle;
    ref_ptr<RelativisticInteraction> relInteraction;
    
    bool haveSecondaries;
    double limit;
    //double thinning;
    mutable std::string interactionTag = "NuANuI";
    
    double neutrinoFieldMass;
    
    std::vector<std::vector<double>> tabEnergy; // 4 columns table depending on the neutrino flavour and number, in initRate they should be built properly
    std::vector<std::vector<double>> tabRate;
    
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    
    std::unordered_map<std::string, int> ratesDictionary;
    mutable std::vector<std::vector<double>> channelProbability; // to be sync with ratesDictionary
    
public:
    
    NeutrinoAntineutrinoInteraction(ref_ptr<NeutrinoField>, ref_ptr<Channels> channels, bool haveSecondaries, double limit, ref_ptr<NeutrinoMixing> neutrinoMixing); // double thinning = 0,
    
    // set the target neutrino field
    void setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField);
    
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
    
    /** Activate/deactivate the neutrino interaction channels
     * @param channels flags
     */
    void setChannels(ref_ptr<Channels> channels);
    void setChannelsBundle(ref_ptr<Channels> channels, std::string fname);
    void setRelativisticInteraction(double m1, double m2, double E, double s);
    
    /** set a custom interaction tag to trace back this interaction
     * @param tag string that will be added to the candidate and output
     */
    void setInteractionTag(std::string tag) const;
    std::string getInteractionTag() const;
    
    void process(crpropa::Candidate *candidate) const;
    void performInteraction(Candidate *candidate, double mass) const;
    
};

} // end namespace crpropa
#endif
