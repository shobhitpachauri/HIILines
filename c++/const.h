#ifndef CONST_INCLUDED
#define CONST_INCLUDED

#include <iostream>
#include <cmath>
#include <gsl/gsl_const_cgsm.h>

const double c = GSL_CONST_CGSM_SPEED_OF_LIGHT;   // km/s
const double h = GSL_CONST_CGSM_PLANCKS_CONSTANT_H;  // J s
const double kb = GSL_CONST_CGSM_BOLTZMANN;    // J/K

const double ergps2Lsun = 1.0 / (3.826e33); // ergs/s to L_sun
const double Jps2Lsun = 1.0 / (3.826e26);    // J/s to L_sun
const double ergs2Jy = 1.0e23;               // erg/s/cm^2/Hz to Jy
const double eV2J = 1.60218E-19;
const double J2erg = 1.0e7;
const double Mpc2cm = 3.086E24;
const double Mpc2km = 3.086E19;
const double kpc2cm = 3.086E21;
const double pc2cm = 3.086E18;
const double Lsun2Jps = 3.826e26;
const double Gyr2s = 1.0e9 * 365 * 24 * 3600;
const double Msun2kg = 1.989E30;
const double mH2kg = 1.6735575e-27; // H mass in unit of kg
const double pi = 3.14159265358979323846;

#endif // CONST_INCLUDED
