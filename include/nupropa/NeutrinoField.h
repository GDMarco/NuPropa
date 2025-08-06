#ifndef NUPROPA_NEUTRINOFIELD_H
#define NUPROPA_NEUTRINOFIELD_H

#include <crpropa/Common.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>

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
        this->mass = 0;
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
    
    int getMass() const {
        return this->mass;
    }
    
protected:

    std::string fieldName;
    bool isRedshiftDependent;
    double mass;
    
};

/**
 @class BlackbodyNeutrinoField
 @brief Neutrino field decorator for black body neutrino fields.
 */
class BlackbodyNeutrinoField: public NeutrinoField {
public:
	BlackbodyNeutrinoField(const std::string fieldName, const double blackbodyTemperature, double mass);
	double getNeutrinoDensity(double eNeutrino, double z = 0.) const;
	double getMinimumNeutrinoEnergy(double z) const;
	double getMaximumNeutrinoEnergy(double z) const;
	void setQuantile(double q);

protected:
	double blackbodyTemperature;
	double quantile;
};

/**
 @class CnuBm1
 @brief Cosmic  neutrino background field with masses m1

 Source info:
 This field is an isotropic blackbody electronic neutrino field with temperature T = 1.93 K
 */
class CnuBm1: public BlackbodyNeutrinoField {
public:
	CnuBm1() : BlackbodyNeutrinoField("CnuBm1", 1.93, 0) {}
};

/**
 @class CnuBm
 @brief Cosmic neutrino background field m2

 Source info:
 This field is an isotropic blackbody electronic neutrino field with temperature T = 1.93 K
 */
class CnuBm2: public BlackbodyNeutrinoField {
public:
    CnuBm2() : BlackbodyNeutrinoField("CnuBm2", 1.93, 8.3e-3 * eV) {}
};

/**
 @class CnuBm3
 @brief Cosmic neutrino background field with mass m3

 Source info:
 This field is an isotropic blackbody neutrino field with temperature T = 1.93 K
 */
class CnuBm3: public BlackbodyNeutrinoField {
public:
    CnuBm3() : BlackbodyNeutrinoField("CnuBm3", 1.93, 50e-3 * eV) {}
};

} // namespace nupropa

#endif // NUPROPA_NEUTRINOFIELD_H
