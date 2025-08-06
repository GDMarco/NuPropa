#include "nupropa/NeutrinoPhotonInteraction.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/ParticleData.h"
#include "nupropa/NeutrinoMixing.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Module.h>
#include <crpropa/Candidate.h>
#include <crpropa/PhotonBackground.h>

#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <filesystem>
#include <unordered_map>

namespace nupropa {

using namespace crpropa;

NeutrinoPhotonInteraction::NeutrinoPhotonInteraction(ref_ptr<PhotonField> photonField, ref_ptr<NeutrinoMixing> neutrinoMixing, bool haveSecondaries,  double limit) { //double thinning,
    
    setPhotonField(photonField);
    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    //setThinning(thinning);
    
    setNeutrinoMixing(neutrinoMixing);
    
}

void NeutrinoPhotonInteraction::setPhotonField(ref_ptr<PhotonField> photonField) {
    
    this->photonField = photonField;
    std::string fname = photonField->getFieldName();
    setDescription("NeutrinoPhotonInteraction::Module" + fname);
    
    std::string pathPh = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoPhotonInteraction/";
    
    initRate(pathPh);
    initCumulativeRate(pathPh);
    
}

// set the neutrino mixing parameters
void NeutrinoPhotonInteraction::setNeutrinoMixing(ref_ptr<NeutrinoMixing> neutrinoMixing) {
    this->neutrinoMixing = neutrinoMixing;
}

void NeutrinoPhotonInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoPhotonInteraction::setLimit(double limit) {
    this->limit = limit;
}

/**
void NeutrinoPhotonInteraction::setThinning(double thinning) {
    this->thinning = thinning;
}
*/

void NeutrinoPhotonInteraction::loadRateFile(const std::string& fileName) {
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

void NeutrinoPhotonInteraction::initRate(std::string filePath) {
    
    tabEnergy.clear();
    tabRate.clear();
    
    std::vector<std::string> flavours = {"Electron", "Muon", "Tauon"};
    std::vector<std::string> masses = {"m1", "m2", "m3"};
    
    int i = 0;
    std::unordered_map<std::string, int> ratesDict;
    
    // 3 masses for each redshift and each flavour! m1, m2, m1, ...
    /**
     Electron .. _m1_z0     0
     Muon ..  m1_z0     1
     Tauon .. m2_z0     2
    ...
    pathPh = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoPhotonInteraction/"; "NeutrinoElectronPhotonInteraction"  " rate_CMB_m1_z0.0.txt"
    */
    for (const auto& f : flavours) {
        for (const auto& a : masses) {
            
            loadRateFile(filePath + "Neutrino" + f + "PhotonInteraction/" +  a + ".txt");
            ratesDictionary[f + "_" + this->photonField->getFieldName() + "_" + a] = i;
            i = i + 1;
            
        }
    }
    
    this->ratesDictionary = ratesDict;
    
}

void NeutrinoPhotonInteraction::loadCumulativeRateFile(const std::string& fileName) {
    
    std::ifstream infile(fileName.c_str());
    
    if (!infile.good())
        throw std::runtime_error("NeutrinoPhotonInteraction: could not open file" + fileName);
    
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

void NeutrinoPhotonInteraction::initCumulativeRate(std::string filePath) {
    
    tabE.clear();
    tabs.clear();
    tabCDF.clear();
    
    std::vector<std::string> flavours = {"Electron", "Muon", "Tauon"};
    std::vector<std::string> masses = {"m1", "m2", "m3"};
    
    // 3 masses for each redshift and each flavour! m1, m2, m1, ...
    /**
     Electron .. _m1_z0     0
     Muon ..  m1_z0     1
     Tauon .. m2_z0     2
    ...
    pathPh = "/Applications/CRPropa/NuNuInteractionv1/CRPropa3-data-zDep/dataOff/NeutrinoInteractions/NeutrinoPhotonInteraction/"; "NeutrinoElectronPhotonInteraction"  " cdf_CMB_m1_z0.0.txt"
    */
    for (const auto& f : flavours) {
        for (const auto& a : masses) {
                
            loadCumulativeRateFile(filePath + "Neutrino" + f + "PhotonInteraction/" +  a + ".txt");
        
        }
    }
}

void NeutrinoPhotonInteraction::setRelativisticInteraction(double m1, double E, double s) const {
    this->relInteraction = new RelativisticInteraction(m1, E, s);
}

int NeutrinoPhotonInteraction::interactionIndex(int ID, double mass) const {
    
    int indexMass = this->neutrinoMixing->massToIndexMass(mass / eV) + 1;
    std::string massComb = this->photonField->getFieldName() + "_m" + std::to_string(indexMass);
    
    std::string alpha;
    if (abs(ID) == 12) {
        alpha = "Electron";
    } else if (abs(ID) == 12) {
        alpha = "Muon";
    } else {
        alpha = "Tau";
    }
    
    auto it = this->ratesDictionary.find(alpha + "_" + massComb);

    int index;
    if (it != this->ratesDictionary.end()) {
        index = it->second;
    } else {
        throw std::runtime_error("Index not found in the dictionary rate!");
    }
    
    return index;
}

class NeutrinoPhotonSecondariesDistribution {
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
        double dcosth;
        double dls;
    
    public:
        
    NeutrinoPhotonSecondariesDistribution(std::string variable, int idChannel, double sThr) {
        
        Ns = 1000;
        Nrer = 1000;
        
        s_min = sThr;
        s_max = 1e28 * eV * eV;
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
        
        if (variable != "costh13_com")
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
            data_i[0] = getDifferentialXS(s, variable, costh_min, idChannel, seedDiffXS) * (dcosth * i);  // to check * (dcosth[i] - 1)!!!
            
            
            for (size_t j = 1; j < Nrer; j++) {
                double costh13_com = costh_min + (j + 0.5) * dcosth;
                double dcosth13_com = (j + 1) * dcosth - j * dcosth; // since we are on a linear scale, dcosth13_com == dcosth
                
                data_i[j] = getDifferentialXS(s, variable, costh13_com, idChannel, seedDiffXS) * dcosth13_com;
                data_i[j] += data_i[j - 1];
            }
            data[i] = data_i;
        }
    }
    
    double getDifferentialXS(double s, std::string variable, double variableValue, int idChannel, int seedDiffXS) {
        
        std::string partonicPath = "/Applications/CRPropa/NuPropaLap/PartonicCalculation/sigmaNu_interface/"; // to change with NUPROPA path
        std::string interfacePath = partonicPath + "bin/";

        std::ostringstream cmd;
        cmd << interfacePath << "Main_Interface.exe"
            << " -c " << idChannel
            << " -s " << seedDiffXS
            << " -E " << std::sqrt(s / GeV / GeV); // Ecms has to be given in GeV
        
        int result = std::system(cmd.str().c_str());
        if (result != 0) {
            std::ostringstream oss;
            oss << "Error: Failed to run Main_Interface.exe, with exit code: " << result;
            throw std::runtime_error(oss.str());
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
            std::ostringstream oss;
            oss << "Error: could not open the differentialXS file. The filename is: " << filename;
            throw std::runtime_error(oss.str());
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

    // draw random costh13_com from the differential cross section distribution
    double sample(double s) {
        size_t idx = std::lower_bound(s_values.begin(), s_values.end(), s) - s_values.begin();
        std::vector<double> s0 = data[idx];
        
        Random &random = Random::instance();
        size_t j = random.randBin(s0); // draw random bin (lower bin boundary returned)
        double costh13_comRand = costh_min + j * dcosth;
        
        return costh13_comRand;
    }
    
};

int NeutrinoPhotonInteraction::fromIDtoChannel(int ID) const {
    if (ID == 12 || ID == 14 || ID == 16)
        return (ID + 44) / 2;   // neutrinos: 12→28, 14→29, 16→30
    if (ID == -12 || ID == -14 || ID == -16)
        return (50 - ID) / 2;   // antineutrinos: -12→31, -14→32, -16→33
    return -1; // invalid ID
}

void NeutrinoPhotonInteraction::performInteraction(Candidate *candidate, int index, double mass) const {
    
    candidate->setActive(false);
    
    if (not haveSecondaries)
        return;
    
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    int ID = candidate->current.getId();
    double w = 1.; // no weights so far!
    
    std::vector<double> vecE = this->tabE[index];
    std::vector<double> vecs = this->tabs[index];
    std::vector<std::vector<double>> vecCDF = this->tabCDF[index];
    
    // check if in tabulated energy range
    if (E < vecE.front() or (E > vecE.back()))
        return;
    
    int leptonID;
    int WbosonID;
    
    int channelID = fromIDtoChannel(ID);
    
    if (ID > 0) {
        leptonID = ID - 1;
        WbosonID = 24; // W+
    } else {
        leptonID = ID + 1;
        WbosonID = -24; // W-
    }
    
    ParticleData particle;
    double leptonMass = particle.getParticleMass(leptonID); // in kg
    double WbosonMass = particle.getParticleMass(WbosonID); // in kg
    
    double sThr = (leptonMass * leptonMass + WbosonMass * WbosonMass) * c_squared * c_squared;
    
    // sample the value of s
    Random &random = Random::instance();
    size_t i = closestIndex(E, vecE);  // find closest tabulation point
    size_t j = random.randBin(vecCDF[i]);
    
    double lo = std::max(sThr, vecs[j-1]); // first s-tabulation point below min(s_kin);
    double hi = vecs[j];
    double s = lo + random.rand() * (hi - lo); // should I add the neutrino masses? since it is the s_kin!!
    
    // sample the cosine of theta13_com
    std::string variable = "costh13_com";
    static NeutrinoPhotonSecondariesDistribution distribution(variable, channelID, sThr);
    double costh13_com = distribution.sample(s);
    
    setRelativisticInteraction(mass, E, s);
    
    // energies of the secondary particles
    std::vector<double> energies = this->relInteraction->getProductEnergiesLab(s, costh13_com, leptonMass, WbosonMass);
  
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    
    //if (haveSecondaries)?, add thinning
    candidate->addSecondary(leptonID, energies[0] / (1 + z), pos, w, interactionTag);
    // W production, to check if on shell?
    candidate->addSecondary(WbosonID, energies[1] / (1 + z), pos, w, interactionTag);
    
}

void NeutrinoPhotonInteraction::process(Candidate *candidate) const
{
    
    // scale the electron energy instead of background photons
    double z = candidate->getRedshift();
    double E = (1 + z) * candidate->current.getEnergy();
    double ID = candidate->current.getId();
    
    if (!(abs(ID) == 12 || abs(ID) == 14 || abs(ID) == 16))
        return;
   
     double mass = this->neutrinoMixing->fromFlavourToMass(ID) * eV; // returned in eV from the function
    
    int index = interactionIndex(ID, mass);
    
    std::vector<double> vecEnergy = this->tabEnergy[index];
    std::vector<double> vecRate = this->tabRate[index];
     
    // check if in tabulated energy range
    if (E < vecEnergy.front() or (E > vecEnergy.back()))
        return;

    // interaction rate
    double rate = interpolate(E, vecEnergy, vecRate);
    rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

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

void NeutrinoPhotonInteraction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string NeutrinoPhotonInteraction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace nupropa
