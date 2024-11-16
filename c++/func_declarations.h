// Filename: func_declarations.h


#ifndef FUNC_DECLARATIONS_H
#define FUNC_DECLARATIONS_H

#include <vector>
#include <iostream>
#include <cmath>
#include <gsl/gsl_const_cgsm.h>


// Hydrogen and helium recombination coefficients
double alphaB_HI(double T4);
double alphaB_HeI(double T4);
double alphaB_Halpha(double T4);
double alphaB_Hbeta(double T4);
double alpha1_HeI(double T4);
double alphaA_HeI(double T4);

// Oxygen recombination coefficients
double alpha_RR_OIII(double T4);
double alpha_RR_OII(double T4);
double alpha_DR_OIII(double T4);
double alpha_DR_OII(double T4);
double alphaB_OIII(double T4);
double alphaB_OII(double T4);

double delta_OII = 1.05e-9;  // cm^3 s^-1
double delta_OI = 1.04e-9;   // cm^3 s^-1

// Charge transfer rate coefficients for Oxygen
double k0_OI_ct(double T4);
double k1_OI_ct(double T4);
double k2_OI_ct(double T4);
double k0r_OI_ct(double T4);

// Photoionization cross sections
std::vector<double> sigma_HI(const std::vector<double>& nu);
std::vector<double> sigma_HeI(const std::vector<double>& nu);
std::vector<double> sigma_HeII(const std::vector<double>& nu);
std::vector<double> sigma_OI(const std::vector<double>& nu);
std::vector<double> sigma_OII(const std::vector<double>& nu);
std::vector<double> sigma_NI(const std::vector<double>& nu);
std::vector<double> sigma_NII(const std::vector<double>& nu);
std::vector<double> sigma_NIII(const std::vector<double>& nu);

// Collisional (de-)excitation coefficients for HeI
double k31_HeI(double T4);
double k32_HeI(double T4);
double k30_HeI(double T4);
double p(double ne, double T4);

// Collisional (de-)excitation coefficients for OII
double Omega10_OII(double T4);
double k10_OII(double T4);
double k01_OII(double T4);
double Omega20_OII(double T4);
double k20_OII(double T4);
double k02_OII(double T4);
double Omega21_OII(double T4);
double k21_OII(double T4);
double k12_OII(double T4);
double Omega30_OII(double T4);
double k30_OII(double T4);
double k03_OII(double T4);
double Omega31_OII(double T4);
double k31_OII(double T4);
double k13_OII(double T4);
double Omega32_OII(double T4);
double k32_OII(double T4);
double k23_OII(double T4);
double Omega40_OII(double T4);
double k40_OII(double T4);
double k04_OII(double T4);
double Omega41_OII(double T4);
double k41_OII(double T4);
double k14_OII(double T4);
double Omega42_OII(double T4);
double k42_OII(double T4);
double k24_OII(double T4);
double Omega43_OII(double T4);
double k43_OII(double T4);
double k34_OII(double T4);
double R01_OII(double ne, double T4);
double R02_OII(double ne, double T4);
double R03_OII(double ne, double T4);
double R04_OII(double ne, double T4);
double R10_OII(double ne, double T4);
double R12_OII(double ne, double T4);
double R13_OII(double ne, double T4);
double R14_OII(double ne, double T4);
double R20_OII(double ne, double T4);
double R21_OII(double ne, double T4);
double R23_OII(double ne, double T4);
double R24_OII(double ne, double T4);
double R30_OII(double ne, double T4);
double R31_OII(double ne, double T4);
double R32_OII(double ne, double T4);
double R34_OII(double ne, double T4);
double R40_OII(double ne, double T4);
double R41_OII(double ne, double T4);
double R42_OII(double ne, double T4);
double R43_OII(double ne, double T4);

// … Other Omega functions for OIII …
// … Other k functions for OIII…
// … Rates R for OIII…

// NII collisional (de-)excitation coefficients
double Omega10_NII(double T4);
double k10_NII(double T4);
double k01_NII(double T4);
double Omega21_NII(double T4);
double Omega31_NII(double T4);
double Omega32_NII(double T4);
double k21_NII(double T4);
double k31_NII(double T4);
double k32_NII(double T4);
double R21_NII(double ne, double T4);
double R31_NII(double ne, double T4);
double R32_NII(double ne, double T4);
// … Other Omega functions for NII …
// … Other k functions for NII…
// … Rates R for NII…

#endif // FUNC_CONST_ION_H