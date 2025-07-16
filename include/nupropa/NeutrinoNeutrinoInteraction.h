#ifndef NUPROPA_NEUTRINONEUTRINOINTERACTION_H
#define NUPROPA_NEUTRINONEUTRINOINTERACTION_H

#include <crpropa/Module.h>
#include "nupropa/NeutrinoBackground.h"
#include "nupropa/RelativisticInteraction.h"

#include <string>
#include <unordered_map>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace nupropa {

using namespace crpropa;

class NeutrinoNeutrinoInteraction : public Module
{
private:
    
    ref_ptr<NeutrinoField> neutrinoField;
    ref_ptr<NeutrinoMixing> neutrinoMixing;
    bool haveSecondaries;
    double limit;
    std::string interactionTag = "NuNuI";
    
    std::unordered_map<std::string, int> ratesDictionary;
    
    double neutrinoFieldMass;
    std::vector<std::vector<double>> tabEnergy;
    std::vector<std::vector<double>> tabRate;
    
    std::vector<std::vector<double>> tabE;
    std::vector<std::vector<double>> tabs;
    std::vector<std::vector<std::vector<double>>> tabCDF;

    ref_ptr<RelativisticInteraction> relInteraction;
    
public:
    
    NeutrinoNeutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, bool haveSecondaries = false, double limit = 0.1, ref_ptr<NeutrinoMixing> neutrinoMixing);
    
    // set the target neutrino field
    void setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField);
    
    // set the neutrino mixing parameters
    void setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing);
    
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
    
    void setRelativisticInteraction(double m1, double m2, double E, double s);
    
    double getDifferentialXS(double s, std::string variable, double variableValue, int idChannel, int seedDiffXS);
    
    void loadRateFile(const std::string& fileName);
    void loadCumulativeRateFile(const std::string& fileName);
    
    void initRate(std::string fileNuNu, std::string fileNuiNuj);
    void initCumulativeRate(std::string fileNuNu, std::string fileNuiNuj);
    
    double findClosestRedshift(double z, const std::vector<double> &redshifts) const;
    int interactionIndex(int ID, int IDbkg, double mass, double z) const;
    
    void process(Candidate *candidate) const;
    void performInteraction(Candidate *candidate, int index, double mass) const;
    
};

} // end namespace nupropa

#endif
