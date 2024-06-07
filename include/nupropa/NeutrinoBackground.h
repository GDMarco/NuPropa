#ifndef NUPROPA_NEUTRINOBACKGROUND_H
#define NUPROPA_NEUTRINOBACKGROUND_H

#include <crpropa/Common.h>
#include <crpropa/Referenced.h>

#include <vector>
#include <string>

namespace nupropa {

using namespace crpropa;
/**
 * \addtogroup NeutrinoFields
 * @{
 */

/**
 @class NeutrinoField
 @brief Abstract base class for neutrino fields.
 */
class NeutrinoField: public Referenced {
public:
    NeutrinoField() {
        this->fieldName = "AbstractNeutrinoField";
        this->isRedshiftDependent = false;
        this->particleID = 12;
    }
    
    /**
     returns comoving neutrino density [1/m^3].
     multiply with (1+z^3) for physical number density.
     @param eNeutrino       neutrino energy [J]
     @param z			redshift (if redshift dependent, default = 0.)
     */
    virtual double getNeutrinoDensity(double eNeutrino, double z = 0.) const = 0;
    virtual double getMinimumNeutrinoEnergy(double z) const = 0;
    virtual double getMaximumNeutrinoEnergy(double z) const = 0;
    virtual std::string getFieldName() const {
        return this->fieldName;
    }
    
    /**
     returns overall comoving scaling factor
     (cf. CRPropa3-data/calc_scaling.py)
     @param z		redshift
     */
    virtual double getRedshiftScaling(double z) const {
        return 1.;
    };
    
    bool hasRedshiftDependence() const {
        return this->isRedshiftDependent;
    }
    
    void setFieldName(std::string fieldName) {
        this->fieldName = fieldName;
    }
    
    int getParticleID() const {     // no needs for being virtual
        return this->particleID;
    }
    
protected:
    
    std::string fieldName;
    bool isRedshiftDependent;
    int particleID;
    
};

/**
 @class TabularNeutrinoField
 @brief Neutrino field decorator for tabulated neutrino fields.

 This class reads neutrino field data from files;
 The first file must be a list of neutrino energies [J], named fieldName_neutrinoEnergy.txt
 The second file must be a list of comoving photon field densities [1/m^3], named fieldName_neutrinoDensity.txt
 Optionally, a third file contains redshifts, named fieldName_redshift.txt

class TabularNeutrinoField: public NeutrinoField {
public:
	TabularNeutrinoField(const std::string fieldName, const bool isRedshiftDependent = true, const int particleID);

	double getNeutrinoDensity(double ePhoton, double z = 0.) const;
	double getRedshiftScaling(double z) const;
	double getMinimumNeutrinoEnergy(double z) const;
	double getMaximumNeutrinoEnergy(double z) const;

protected:
	void readNeutrinoEnergy(std::string filePath);
	void readNeutrinoDensity(std::string filePath);
	void readRedshift(std::string filePath);
	void initRedshiftScaling();
	void checkInputData() const;

	std::vector<double> photonEnergies;
	std::vector<double> photonDensity;
	std::vector<double> redshifts;
	std::vector<double> redshiftScalings;
};
 */
 
/**
 @class IRB_Kneiske04
 @brief Extragalactic background light model from Kneiske et al. 2004

 Source info:
 DOI:10.1051/0004-6361:20031542,
 https://www.aanda.org/articles/aa/pdf/2004/03/aa3848.pdf, figure 1 ("Best-fit" model)

class IRB_Kneiske04: public TabularPhotonField {
public:
	IRB_Kneiske04() : TabularPhotonField("IRB_Kneiske04", true) {}
};
*/

/**
 @class BlackbodyNeutrinoField
 @brief Neutrino field decorator for black body neutrino fields.
 */
class BlackbodyNeutrinoField: public NeutrinoField {
public:
	BlackbodyNeutrinoField(const std::string fieldName, const double blackbodyTemperature, int particleID);
	double getNeutrinoDensity(double eNeutrino, double z = 0.) const;
	double getMinimumNeutrinoEnergy(double z) const;
	double getMaximumNeutrinoEnergy(double z) const;
	void setQuantile(double q);

protected:
	double blackbodyTemperature;
	double quantile;
};

/**
 @class CnuBe
 @brief Cosmic electronic neutrino background field

 Source info:
 This field is an isotropic blackbody electronic neutrino field with temperature T = 1.93 K
 */
class CnuBe: public BlackbodyNeutrinoField {
public:
	CnuBe() : BlackbodyNeutrinoField("CnuBe", 1.93, 12) {}
};

/**
 @class CnuBm
 @brief Cosmic muonic neutrino background field

 Source info:
 This field is an isotropic blackbody electronic neutrino field with temperature T = 1.93 K
 */
class CnuBm: public BlackbodyNeutrinoField {
public:
    CnuBm() : BlackbodyNeutrinoField("CnuBm", 1.93, 14) {}
};

/**
 @class CnuBt
 @brief Cosmic tauonic neutrino background field

 Source info:
 This field is an isotropic blackbody tauonic neutrino field with temperature T = 1.93 K
 */
class CnuBt: public BlackbodyNeutrinoField {
public:
    CnuBt() : BlackbodyNeutrinoField("CnuBt", 1.93, 16) {}
};

/**
 @class CnuxBe
 @brief Cosmic electronic antineutrino background field

 Source info:
 This field is an isotropic blackbody electronic antineutrino field with temperature T = 1.93 K
 */
class CnuxBe: public BlackbodyNeutrinoField {
public:
    CnuxBe() : BlackbodyNeutrinoField("CnuxBe", 1.93, -12) {}
};

/**
 @class CnuxBm
 @brief Cosmic muonic antineutrino background field

 Source info:
 This field is an isotropic blackbody electronic antineutrino field with temperature T = 1.93 K
 */
class CnuxBm: public BlackbodyNeutrinoField {
public:
    CnuxBm() : BlackbodyNeutrinoField("CnuxBm", 1.93, -14) {}
};

/**
 @class CnuxBt
 @brief Cosmic tauonic antineutrino background field

 Source info:
 This field is an isotropic blackbody tauonic antineutrino field with temperature T = 1.93 K
 */
class CnuxBt: public BlackbodyNeutrinoField {
public:
    CnuxBt() : BlackbodyNeutrinoField("CnuxBt", 1.93, -16) {}
};

} // namespace nupropa

#endif // NUPROPA_NEUTRINOBACKGROUND_H
