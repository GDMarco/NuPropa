#ifndef NUPROPA_NEUTRINONEUTRINOINTERACTION_H
#define NUPROPA_NEUTRINONEUTRINOINTERACTION_H

#include <crpropa/Module.h>
#include <crpropa/Random.h>
#include "nupropa/NeutrinoField.h"
#include "nupropa/NeutrinoMixing.h"
#include "nupropa/RelativisticInteraction.h"

#include <iomanip>
#include <string>
#include <unordered_map>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <sys/stat.h>

namespace nupropa {

using namespace crpropa;

// the structure is inspired by EMInverseCompton class etc.
class NeutrinoNeutrinoSecondariesDistribution {
    private:
        // channelId → [sIndex][cosIndex] cumulative table
        std::unordered_map<int, std::vector<std::vector<double>>> dataPerChannel;
    
        std::vector<double> s_values;
        std::vector<double> costh13_com_values;
        size_t Ns = 1000;
        size_t Nrer = 1000;

        double s_min = 1e-100;
        double s_max = 1e28 * eV * eV;
        double costh_min = -1.0;
        double costh_max = 1.0;
        double dcosth;
        double dls;

        std::string variable;
    
    public:
    
    explicit NeutrinoNeutrinoSecondariesDistribution(const std::string& variable_) : variable(variable_) {
        
        if (variable != "dsigdcosth")
            throw std::runtime_error("Only dsigdcosth is supported");
        
    }

    void buildChannel(int idChannel) {
        
        if (dataPerChannel.count(idChannel))
            return; // already built
        
        dls = (log(s_max) - log(s_min)) / Ns;
        dcosth = (costh_max - costh_min) / Nrer;
        
        // s grid
        s_values = std::vector<double>(1001);
        for (size_t i = 0; i <= Ns; ++i)
            s_values[i] = s_min * std::exp(i * dls);
        
        // costh grid
        costh13_com_values = std::vector<double>(1001);
        for (size_t i = 0; i <= Nrer; ++i)
            costh13_com_values[i] = costh_min + i * dcosth;
        
        std::vector<std::vector<double>> data(Ns, std::vector<double>(Nrer));
        std::vector<double> data_i(Nrer);
        
        int seedDiffXS = 1;
        
        for (size_t i = 0; i < Ns; ++i) {
            
            double s = s_min * std::exp((i + 0.5) * dls);
            
            for (size_t j = 0; j < Nrer; ++j) {
                
                double costh13_com = costh_min + (j + 0.5) * dcosth;
                
                double val =
                getDifferentialXS(s, variable, costh13_com, idChannel, seedDiffXS)
                * dcosth;
                
                if (j == 0) {
                    data_i[j] = val;
                } else {
                    data_i[j] = val + data_i[j - 1];
                }
            }
            
            data[i] = data_i;
        }
        
        dataPerChannel.emplace(idChannel, std::move(data));
        
    }
    
        
    static std::string resolvePartonicPath() {
        const char *env = std::getenv("NUPROPA_PARTONIC_PATH");
        std::string path = (env && *env)
            ? std::string(env)
            : "/sdf/home/g/gaetano/CRPropa/NuPropa/PartonicCalculation/sigmaNu_interface/";
        if (!path.empty() && path.back() != '/') {
            path += "/";
        }
        return path;
    }

    static bool fileExists(const std::string& path) {
        struct stat st;
        return ::stat(path.c_str(), &st) == 0;
    }

    double getDifferentialXS(double s, std::string variable, double variableValue, int idChannel, int seedDiffXS) {
        
        std::string partonicPath = resolvePartonicPath(); // can be overridden by NUPROPA_PARTONIC_PATH
        std::string interfacePath = partonicPath + "bin/";
        if (std::isnan(s) || std::isinf(s)) {
            std::cerr << "s became invalid at s=" << s << ", varValue: "<< variableValue << std::endl;
            
        }
        
        std::ostringstream ss;
        ss << std::scientific << std::setprecision(5) << std::sqrt(s / GeV / GeV);
        std::string Ecms = ss.str();
        
        double Ecms_val = std::sqrt(s / GeV / GeV);
        
        if (std::isnan(Ecms_val) || !std::isfinite(Ecms_val)) { // to format in CRPropa style!
            std::ostringstream err;
            err << "\n[FATAL] Ecms is NaN or infinite\n"
            << "--------------------------------\n"
            << "s               = " << s << "\n"
            << "GeV             = " << GeV << "\n"
            << "idChannel       = " << idChannel << "\n"
            << "variable        = " << variable << "\n"
            << "variableValue   = " << variableValue << "\n"
            << "seedDiffXS      = " << seedDiffXS << "\n"
            << "s/(GeV*GeV)     = " << s/(GeV*GeV) << "\n";
            
            throw std::runtime_error(err.str());
        }
        
        std::string filePath = partonicPath + "dataDifferentialXS/channel" + std::to_string(idChannel) + "/";
        std::string filename = filePath + variable + "_channel" + std::to_string(idChannel) +
        "_Ecms" + Ecms + "_s" + std::to_string(seedDiffXS) + ".txt";
        
        
        
        const std::string interfaceExe = interfacePath + "Main_Interface.exe";
        if (!fileExists(interfaceExe)) {
            std::ostringstream oss;
            oss << "Error: Main_Interface.exe not found at " << interfaceExe
                << ". Set NUPROPA_PARTONIC_PATH to the sigmaNu_interface directory or build the interface.";
            throw std::runtime_error(oss.str());
        }

        if (!std::ifstream(filename)) {
            std::ostringstream cmd;
            cmd << interfaceExe
            << " -c " << idChannel
            << " -s " << seedDiffXS
            << " -E " << std::sqrt(s / GeV / GeV); // Ecms in GeV
            
            int result = std::system(cmd.str().c_str());
            if (result != 0) {
                std::ostringstream oss;
                oss << "Error: Failed to run Main_Interface.exe, with exit code: " << result;
                throw std::runtime_error(oss.str());
            }
        }
        
        std::vector<double> variableScan;
        std::vector<double> differentialXS;
        
        std::ifstream infile(filename);
        
        if(!infile) {
            std::ostringstream oss;
            oss << "Error: could not open the differentialXS file. The filename is: " << filename;
            throw std::runtime_error(oss.str());
        } else {
            std::string line;
            double a, b;
            
            while (std::getline(infile, line)) {
                // skip empty lines or lines starting with '#'
                if (line.empty() || line[0] == '#') continue;
                
                std::istringstream iss(line);
                if (iss >> a >> b) {
                    variableScan.push_back(a);
                    differentialXS.push_back(b);
                } else {
                    // if parsing failed on a non-comment line, optionally throw or skip
                    std::ostringstream oss;
                    oss << "Error: could not parse line: \"" << line << "\" in file " << filename;
                    throw std::runtime_error(oss.str());
                }
            }
            infile.close();
        }
        
        double diffXS = interpolate(variableValue, variableScan, differentialXS);
        return diffXS;
        
    }
    
    double sample(int idChannel, double s, double mass, double massField) {
        
        if (!dataPerChannel.count(idChannel)) {
            buildChannel(idChannel);
        }
        const auto channelIt = dataPerChannel.find(idChannel);
        if (channelIt == dataPerChannel.end()) {
            throw std::runtime_error("Missing secondary distribution for channel " + std::to_string(idChannel));
        }
        const auto& data = channelIt->second;
        
        double sThr = std::max(std::pow(mass + massField, 2), 1e-100);
        
        // to check!
        auto sBegin = std::lower_bound(s_values.begin(), s_values.end(), sThr);
        auto sIt = std::lower_bound(sBegin, s_values.end(), s);
        size_t idx = std::distance(s_values.begin(), sIt);
        
        // check the validity! + apply other protections + apply to the other modules!
        if (idx >= data.size()) {
            idx = data.size() - 1;
        }
        const std::vector<double>& s0 = data[idx];
        
	Random &random = Random::instance(); 
        size_t j = random.randBin(s0);
        return costh_min + j * dcosth;
    }
        
};

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

    mutable ref_ptr<RelativisticInteraction> relInteraction;
    
    mutable NeutrinoNeutrinoSecondariesDistribution secondariesDistribution;
    
public:
    
    NeutrinoNeutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, ref_ptr<NeutrinoMixing> neutrinoMixing, bool haveSecondaries = false, double limit = 0.1);
    
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
    
    void setRelativisticInteraction(double m1, double m2, double E, double s) const;
    
    void loadRateFile(const std::string& fileName);
    void loadCumulativeRateFile(const std::string& fileName);
    
    void buildSecondariesDistributionClass();
    
    void initRate(std::string fileNuNu, std::string fileNuiNuj);
    void initCumulativeRate(std::string fileNuNu, std::string fileNuiNuj);
    
    double findClosestRedshift(double z, const std::vector<double> &redshifts) const;
    int interactionIndex(int ID, int IDbkg, double mass, double z) const;
    
    void process(Candidate *candidate) const;
    void performInteraction(Candidate *candidate, int index, double mass, int IDbkg) const;
    
};

} // end namespace nupropa

#endif
