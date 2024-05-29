#include <crpropa/Module.h>
#include "NeutrinoAntineutrinoInteraction.h"
#include "Channels.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <string>


namespace crpropa {
/// A custom C++ module for Neutrino-Neutrino Interactions in astorphysical scenarios
class NeutrinoAntineutrinoInteraction : public Module
{
private:
    
    ref_ptr<NeutrinoField> neutrinoField;
    int neutrinoFieldID; // = 12; // check the initialization!
    bool haveSecondaries;
    ref_ptr<Channels> channels;
    double limit;
    //double thinning;
    
    std::string interactionTag = "NuAntiNuI";
    
    std::vector<std::vector<double>> tabEnergy; // 4 columns table depending on the neutrino flavour and number, in initRate they should be built properly
    std::vector<std::vector<double>> tabRate;
    
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;
    std::vector<std::vector<int>> tabProductsID;
    
    std::unordered_map<int, std::string> interactionDictionary;
    std::vector<std::vector<double>> channelProbability; // to be sync with interactionDictionary
public:
    /// The parent's constructor need to be called on initialization!
    NeutrinoNeutrinoInteraction(ref_ptr<NeutrinoField>, bool haveSecondaries = false, ref_ptr<Channels> channels, double limit = 0.1); // double thinning = 0,
    
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
    
    int searchChannel(std::string interacting, std::string products);
    std::vector<double> fillTableZeros(std::vector<double> table, int size); // table is the vector and the size is the size I want to obtain
    // fill the vector with zeros from the beginning
    std::vector<double> getTabulatedRateEnergy(int ID, int nuBkgID);
    void computeInteractionProbabilities (std::vector<std::vector<double>> rates);
    
    /** set a custom interaction tag to trace back this interaction
     * @param tag string that will be added to the candidate and output
     */
    void setInteractionTag(std::string tag);
    std::string getInteractionTag() const;
    
    void initRate(std::string fname);
    // void initCumulativeRate(std::string filename);
    
    void process(crpropa::Candidate *candidate) const;
    void performInteraction(Candidate *candidate) const;
    
};

} // end namespace crpropa
