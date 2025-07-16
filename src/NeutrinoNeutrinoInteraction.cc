#include "nupropa/NeutrinoNeutrinoInteraction.h"
#include "nupropa/NeutrinoBackground.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/NeutrinoMixing.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace nupropa {

using namespace crpropa;

NeutrinoNeutrinoInteraction::NeutrinoNeutrinoInteraction(ref_ptr<NeutrinoField> neutrinoField, bool haveSecondaries, double limit, ref_ptr<NeutrinoMixing> neutrinoMixing) : Module() {
    setNeutrinoField(neutrinoField);
    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    setNeutrinoMixing(neutrinoMixing);
}

void NeutrinoNeutrinoInteraction::setNeutrinoField(ref_ptr<NeutrinoField> neutrinoField) {
    this->neutrinoField = neutrinoField;
    this->neutrinoFieldMass = neutrinoField->getMass();
    std::string fname = neutrinoField->getFieldName();
    setDescription("NeutrinoNeutrinoInteraction::Module" + fname);
    setInteractionTag("NuNuInt");
    
    std::string pathNuNu = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoNeutrinoInteraction/NeutrinoNeutrinoElastic/;
    // getDataPath("NeutrinoNeutrinoInteraction/NeutrinoNeutrinoElastic/"
    std::string pathNuiNuj = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoNeutrinoInteraction/NeutrinoiNeutrinojElastic/;
    // getDataPath("NeutrinoNeutrinoInteraction/NeutrinoiNeutrinojElastic/"
    
    initRate(pathNuNu + "rate_" + fname, pathNuiNuj + "rate_" + fname);
    initCumulativeRate(pathNuNu + "cdf_" + fname, pathNuiNuj + "cdf_" + fname);
    
}

void NeutrinoNeutrinoInteraction::setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    this->neutrinoMixing = neutrinoMixing;
}

void NeutrinoNeutrinoInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoNeutrinoInteraction::setLimit(double limit) {
    this->limit = limit;
}

void NeutrinoNeutrinoInteraction::loadRateFile(const std::string& fileName) {
    std::ifstream infile(fileName.c_str());
    if (!infile.good()) {
        throw std::runtime_error("Could not open file: " + fileName);
    }

    std::vector<double> vecEnergy, vecRate;

    while (infile.good()) {
        if (infile.peek() != '#') {
            double a, b;
            infile >> a >> b;
            if (infile) {
                vecEnergy.push_back(std::pow(10, a) * eV);
                vecRate.push_back(b / Mpc);
            }
        }
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    infile.close();

    this->tabEnergy.push_back(vecEnergy);
    this->tabRate.push_back(vecRate);
}

void NeutrinoNeutrinoInteraction::initRate(std::string pathNuNu, std::string pathNuiNuj) {
    
    this->tabEnergy.clear();
    this->tabRate.clear();
    
    std::vector<std::string> masses = {"m1", "m2", "m3"};
    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50}; // same as computed in NuPropa-data
    
    int i = 0;
    std::unordered_map<std::string, int> ratesDict;
    
    // 9 combination of masses! m1_m1, m1_m2, m2_m1, ...
    /**
     m1_m1_z0 NuNu 0
     m1_m1_z0 NuiNuj 1
     m1_m2_z0 NuNu 2
    ...
     */
    for (const auto& a : masses) {
        for (const auto& z : redshifts) {
            
            std::ostringstream out;
            out << std::fixed << std::setprecision(1) << z;
            std::string zDec = out.str();
            
            loadRateFile(pathNuNu + a + "_z" + zDec + ".txt");
            ratesDictionary[this->neutrinoField->getFieldName + "_" + a + "_z" + zDec + "_NuNu"] = i;
            i = i + 1;
            
            loadRateFile(pathNuiNuj + a + "_z" + zDec + ".txt");
            ratesDictionary[this->neutrinoField->getFieldName + "_" + a + "_z" + zDec + "_NuiNuj"] = i;
            i = i + 1;
        }
    }
    this->ratesDictionary = ratesDict;
}

void NeutrinoNeutrinoInteraction::loadCumulativeRateFile(const std::string& fileName) {
    
    std::ifstream infile(filename.c_str());
    
    if (!infile.good())
        throw std::runtime_error("NeutrinoNeutrinoInteraction: could not open file" + filename);
    
    // skip header
    while (infile.peek() == '#')
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    
    // declare the row vectors
    std::vector<double> vecE;
    std::vector<double> vecs;
    std::vector<std::vector<double>> vecCDF;
    
    // read s values in the first line
    double a;
    infile >> a; // skip the first value
    while (infile.good() and (infile.peek() != '\n')) {
        infile >> a;
        vecs.push_back(pow(10, a) * eV * eV);
    }
    
    // read all the following lines: E, cdf values
    while(infile.good()) {
        infile >> a;
        if (!infile)
            break; // end of file
        
        vecE.push_back(pow(10, a) * eV);
        
        std::vector<double> cdf;
        for (int i = 0; i < vecs.size(); i++) {
            infile >> a;
            cdf.push_back(a / Mpc);
        }
        vecCDF.push_back(cdf);
    }
    
    this->tabs.push_back(vecs);
    this->tabE.push_back(vecE);
    this->tabCDF.push_back(vecCDF);

}

void initCumulativeRate(std::string pathNuNu, std::string pathNuiNuj) {

    this->tabE.clear();
    this->tabs.clear();
    this->tabCDF.clear();
    
    std::vector<std::string> masses = {"m1", "m2", "m3"};
    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50}; // same as computed in NuPropa-data
    
    for (const auto& a : masses) {
        for (const auto& z : redshifts) {
            
            std::ostringstream out;
            out << std::fixed << std::setprecision(1) << z;
            std::string zDec = out.str();
            
            loadCumulativeRateFile(pathNuNu + a + "_z" + zDec + ".txt");
            loadCumulativeRateFile(pathNuiNuj + a + "_z" + zDec + ".txt");
        }
    }
}

void NeutrinoNeutrinoInteraction::setRelativisticInteraction(double m1, double m2, double E, double s) {
    this->relInteraction = new relInteraction(m1, m2, E, s);
}

double NeutrinoNeutrinoInteraction::getDifferentialXS(double s, std::string variable, double variableValue, int idChannel, int seedDiffXS) {
    
    std::string partonicPath = "/Applications/CRPropa/NuPropaLap/PartonicCalculation/sigmaNu_interface/"; // to change with NUPROPA path
    std::string interfacePath = partonicPath + "bin/";

    std::ostringstream cmd;
    cmd << interfacePath << "Main_Interface.exe"
        << " -c " << idChannel
        << " -s " << seedDiffXS
        << " -E " << std::sqrt(s / GeV / GeV); // Ecms has to be given in GeV
    
    int result = std::system(cmd.str().c_str());
    if (result != 0) {
        throw std::runtime_error("Error: Failed to run Main_Interface.exe, with exit code: " << result);
    }

    std::ostringstream ss;
    ss << std::scientific << std::setprecision(5) << s / GeV / GeV;
    std::string Ecms2 = ss.str();
    
    std::string filePath = partonicPath + "dataDifferentialXS/channel" + std::to_string(idChannel) + "/";
    std::string filename = filePath + variable + "_channel" + std::to_string(idChannel) + "_EcmsSq" + Ecms2 + "_s" + std::to_string(seedDiffXS) + ".txt";
    
    std::vector<double> variableScan;
    std::vector<double> differentialXS;
    
    std::ifstream infile(filename);
    
    if(!infile) {
        throw std::runtime_error( << "Error: could not open the diffentialXS file. The filename is: " << filename );
    } else {
        
        double a, b;
        while (infile >> a >> b) {
            variableScan.push_back(a);
            differentialXS.push_back(b);
        }
        infile.close();
    }
    
    double diffXS = interpolate(variableValue, variableScan, differentialXS);
    return diffXS;
    
}

// the structure is inspired by EMInverseCompton class etc.
class NeutrinoNeutrinoSecondariesDistribution {
    private:
        std::vector< std::vector<double> > data;
        std::vector<double> s_values;
        std::vector<double> costh13_com_values;
        size_t Ns;
        size_t Nrer;
        double s_min;
        double s_max;
        double costh_min;
        double costh_max;
        double dls;
    
    public:
        
    NeutrinoNeutrinoSecondariesDistribution(std::string variable, int idChannel, double mass) {
        
        Ns = 1000;
        Nrer = 1000;
        
        s_min = mass * mass + this->neutrinoFieldMass * this->neutrinoFieldMass;
        s_max = 1e28 * eV * eV; // value from the computed tables
        dls = (log(s_max) - log(s_min)) / Ns;
        
        costh_min = -1;
        costh_max = 1;
        dcosth = (costh_max - costh_min) / Nrer;
        
        data = std::vector< std::vector<double> >(1000, std::vector<double>(1000));
        std::vector<double> data_i(1000);
        
        // tabulating the s bin borders
        s_values = std::vector<double>(1001);
        for (size_t i = 0; i < Ns + 1; ++i)
            s_values[i] = s_min * exp(i * dls);
        
        if (!variable == "costh13_com")
            throw std::runtime_error("The only available variable to compute the differential cross section is costheta13_com!");
        
        // tabulating the costh13_com bin borders
        costh13_com_values = std::vector<double>(1001);
        for (size_t i = 0; i < Nrer + 1; ++i)
            costh13_com_values[i] = costh_min + i * dcosth;
        
        // to generate randomly the seed to produce the differential cross section
        // Random &random = Random::instance();
        // int seedDiffXS = random.randInt();
        // not affecting the results
        
        int seedDiffXS = 1;
        
        // for each s and costh13_com tabulate the cumulative differential cross section
        for (size_t i = 0; i < Ns; i++) {
            double s = s_min * exp((i + 0.5) * dls);
            
            // cumulative midpoint integration
            data_i[0] = getDifferentialXS(s, variable, costh_min, idChannel, seedDiffXS) * (dcosth[i] - 1);
            // to check * (dcosth[i] - 1)!!!
            
            for (size_t j = 1; j < Nrer; j++) {
                double costh13_com = costh_min + (j + 0.5) * dcosth;
                double dcosth13_com = (j + 1) * dcosth - j * dcosth; // since we are on a linear scale, dcosth13_com == dcosth
                
                data_i[j] = getDifferentialXS(s, variable, cost13_com, idChannel, seedDiffXS) * dcosth13_com;
                data_i[j] += data_i[j - 1];
            }
            data[i] = data_i;
        }
        
    }
    
    // draw random costh13_com from the differential cross section distribution
    double sample(double s) {
        size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
        std::vector<double> s0 data[idx];
        
        Random &random = Random::instance();
        size_t j = random.randBin(s0); // draw random bin (lower bin boundary returned)
        double costh13_comRand = costh_min + j * dcosth;
        
        return costh13_comRand;
    }
        
};

void NeutrinoNeutrinoInteraction::performInteraction(Candidate *candidate, int index, double mass) const {
    
    // check if in tabulated energy range
    if (E < tabE.front() or (E > tabE.back()))
            return;

    double E = candidate->current.getEnergy();
    double ID = candidate->current.getId();
    double w = 1; // no thinning, TBD
    
    // here I just need the index found before in computeRate
    std::vector<double> vecE = this->tabE[index];
    std::vector<double> vecs = this->tabs[index];
    std::vector<std::vector<double>> vecCDF = this->tabCDF[index];
    
    size_t i = closestIndex(E, vecE);  // find closest tabulation point
    size_t j = random.randBin(vecCDF[i]);
    double lo = std::max(mass * mass + this->neutrinoFieldMass * this->neutrinoFieldMass, tabs[j-1]); // first s-tabulation point below min(s_kin);
    
    double hi = vecs[j];
    double s = lo + random.rand() * (hi - lo); // should I add the neutrino masses? since it is the s_kin!!
    
    int idInteraction;
    
    if (ID == this->neutrinoFieldID) {
        if (ID > 0) {
            idInteraction = 1; // in Rhorry's "nu_e nu_e -> nu_e nu_e"
        } else {
            idInteraction = 3; // in Rhorry's "nu_e~ nu_e~ -> nu_e~ nu_e~"
        }
    } else {
        if (ID > 0) {
            idInteraction = 5; // in Rhorry's "nu_e nu_mu -> nu_e nu_mu"
        } else {
            idInteraction = 8; // in Rhorry's "nu_e~ nu_mu~ -> nu_e~ nu_mu~"
        }
    }
    
    // sample the cosine of theta13_com, the only available for now
    std::string variable = "costh13_com";
    static NeutrinoAntineutrinoSecondariesDistribution distribution(variable, idInteraction, mass);
    double costh13_com = distribution.sample(s);
    
    // see if this function wants the neutrino masses in J or kg!
    setRelativisticInteraction(mass / c_squared, this->neutrinoFieldMass / c_squared, E);
    
    // energies of the secondary particles
    std::vector<double> energies = this->relInteraction->getProductEnergiesLab(s, costh13_com, mass / c_squared, this->neutrinoFieldMass / c_squared);
    
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    
    candidate->current.setEnergy(energies[0]);
    
    // add the thinning
    if (haveSecondaries)
        candidate->addSecondary(this->neutrinoFieldID, energies[1], pos, w, interactionTag);
}

double NeutrinoNeutrinoInteraction::findClosestRedshift(double z, const std::vector<double> &redshifts) const {
    auto it = std::lower_bound(redshifts.begin(), redshifts.end(), z);

    if (it == redshifts.begin()) return *it;
    if (it == redshifts.end()) return redshifts.back();

    double upper = *it;
    double lower = *(it - 1);

    return (std::abs(upper - z) < std::abs(lower - z)) ? upper : lower;
}

int NeutrinoNeutrinoInteraction::interactionIndex(int ID, int IDbkg, double mass, double z) const {
    
    int indexMass = this->mixingNeutrino->massToIndexMass(mass / eV) + 1;
    std::string massComb = this->neutrinoField->getFieldName + "_m" + std::to_string(indexMass);
    
    std::string alphaBeta;
    if (ID == this->neutrinoFieldID) {
        alphaBeta = "_NuNu";
    } else {
        alphaBeta = "_NuiNuj";
    }
    
    std::vector<double> redshifts = {0, 2, 5, 8, 11, 15, 20, 25, 30, 40, 50}; // same as computed in NuPropa-data
    
    double zClosest = findClosestRedshift(z, redshifts);
    std::ostringstream out;
    out << std::fixed << std::setprecision(1) << zClosest;
    std::string zDec = out.str();
    
    std::string redshift = "_z" + zDec;

    auto it = this->ratesDictionary.find(massComb + redshift + alphaBeta);

    int index;
    if (it != this->ratesDictionary.end()) {
        index = it->second;
    } else {
        throw std::runtime_error("Index not found in the dictionary rate!");
    }
    
    return index;
}

void NeutrinoNeutrinoInteraction::process(Candidate *candidate) const {

    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy();
    int ID = candidate->current.getId();
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;
    
    int IDbkg = this->neutrinoMixing->fromMassToFlavour(this->neutrinoFieldMass / eV);
    
    if (!(ID * IDbkg > 0))
        return;
    
    double mass = this->neutrinoMixing->fromFlavourToMass(ID) * eV; // returned in eV from the function
    
    std::vector<double> vecEnergy;
    std::vector<double> vecRate;
    
    int index = interactionIndex(ID, IDbkg, mass, z);
    
    vecEnergy = this->tabEnergy[index];
    vecRate = this->tabRate[index];
    
    // check if in tabulated energy range
    if (E < vecEnergy.front() or (E > vecEnergy.back()))
        return;
        
    // interaction rate
    double rate = interpolate(E, vecEnergy, vecRate);
    
    // check for interaction
    Random &random = Random::instance();
    double randDistance = -log(random.rand()) / rate;
    double step = candidate->getCurrentStep();
    if (step < randDistance) {
        candidate->limitNextStep(limit / rate);
        return;
    } else { // after performing interaction neutrino ceases to exist (hence return)
        performInteraction(candidate, index, mass);
        return;
    }
}

void NeutrinoNeutrinoInteraction::setInteractionTag(std::string tag) {
    this->interactionTag = tag;
}

std::string NeutrinoNeutrinoInteraction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace nupropa

