#ifndef NUPROPA_PARTICLEDATA_H
#define NUPROPA_PARTICLEDATA_H

#include <crpropa/Module.h>

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <unordered_map>

namespace nupropa {

using namespace crpropa;

class ParticleData
{
private:
    
    std::unordered_map<int, double> dictionaryIDmass;
    
public:
    
    ParticleData();
    ParticleData(std::unordered_map<int, double> IDmass);
    
    void setParticleIDmass(std::unordered_map<int, double> IDmass);
    
    double getParticleMass(int ID);
    void addNewParticle(int ID, double mass);
    
};

} // end namespace nupropa

#endif // PARTICLEDATA_H
