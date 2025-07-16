#include "nupropa/NeutrinoBackground.h"
#include <crpropa/Units.h>
#include <crpropa/Random.h>

#include "kiss/logger.h"

#include <fstream>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace nupropa {

using namespace crpropa;
 
BlackbodyNeutrinoField::BlackbodyNeutrinoField(std::string fieldName, double blackbodyTemperature, double mass) {
	this->fieldName = fieldName;
	this->blackbodyTemperature = blackbodyTemperature;
    this->mass = mass,
	this->quantile = 0.0001; // tested to be sufficient, only used for extreme values of primary energy or temperature
}

double BlackbodyNeutrinoField::getNeutrinoDensity(double Eneutrino, double z) const {
	return 8 * M_PI * pow_integer<3>(Eneutrino / (h_planck * c_light)) / (std::exp(Eneutrino / (k_boltzmann * this->blackbodyTemperature)) + 1) / 2; // the /2 is because the BB density should be for the nu+nux, SEE better!
}

double BlackbodyNeutrinoField::getMinimumNeutrinoEnergy(double z) const {
	double A;
	int quantile_int = 10000 * quantile;
	switch (quantile_int)
	{
	case 1:	// 0.01 % percentil
		A = 1.093586e-5 * eV / kelvin;
		break;
	case 10:		// 0.1 % percentil
		A = 2.402189e-5 * eV / kelvin;
		break;
	case 100:		// 1 % percentil
		A = 5.417942e-5 * eV / kelvin;
		break;
	default:
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
		break;
	}
	return A * this -> blackbodyTemperature;
}

double BlackbodyNeutrinoField::getMaximumNeutrinoEnergy(double z) const {
	double factor = std::max(1., blackbodyTemperature / 2.73);
	return 0.1 * factor * eV; // T dependent scaling, starting at 0.1 eV as suitable for CMB
}

void BlackbodyNeutrinoField::setQuantile(double q) {
	if(not ((q == 0.0001) or (q == 0.001) or (q == 0.01)))
		throw std::runtime_error("Quantile not understood. Please use 0.01 (1%), 0.001 (0.1%) or 0.0001 (0.01%) \n");
	this -> quantile = q;
}

} // namespace crpropa
