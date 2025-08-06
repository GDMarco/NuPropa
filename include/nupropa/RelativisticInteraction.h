#ifndef NUPROPA_RELATIVISTICINTERACTION_H
#define NUPROPA_RELATIVISTICINTERACTION_H

#include <crpropa/Referenced.h>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>

namespace nupropa {

using namespace crpropa;

/// A custom C++ class for the treatment of the interaction in the center of mass frame.
class RelativisticInteraction : public Referenced {

private:
    
    double eps;
    double beta_com;
    double gamma_com;
    
public:
    
    RelativisticInteraction();
    
    RelativisticInteraction(double m1, double m2, double E, double s);
    
    RelativisticInteraction(double m1, double E, double s);
    
    void setBetaCom(double E, double m1, double m2, double s);
    void setBetaPhotonCom(double E, double m1, double s);
    
    void setGammaCom(double E, double s);
    
    double getBetaCom() {
        return this->beta_com;
    }
    double getGammaCom() {
        return this->gamma_com;
    }
    
    double computeProductsMomentumCom(double s, double m3, double m4);
    
    std::vector<double> getProductEnergiesLab(double s, double costh13_com, double m3, double m4);
    
};
    
} // end namespace nupropa

#endif // RELATIVISTICINTERACTION_H
