#pragma once
#include <cmath>
#include <complex>
#include "tools.hh"

// Master function to evaluate scattering processes, per channel
extern double dsigma_channels(KinematicData &, int );


// Implementation of inclusive result for channels 1,2,6 (corresponding to four neutrino scattering)
extern double sigma_nu_incl(double shat, int chan);

extern double sigma_Wl_incl(double shat, int chan);

extern double sigma_ZZ_incl(double shat);
extern double sigma_ZZ_incl_Rhorry(double shat);
extern double sigma_WW_incl(double shat);
