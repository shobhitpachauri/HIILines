#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "const.h"
#include <gsl/gsl_const_cgsm.h>
#include "func_declarations.h"

// Recombination coefficients
double alphaB_HI(double T4) {
    return 2.59e-13 * pow(T4, (-0.833 - 0.034 * log(T4))); // cm^3 s^-1
}

double alphaB_HeI(double T4) {
    return 2.72e-13 * pow(T4, -0.789); // cm^3 s^-1
}

double alphaB_Halpha(double T4) {
    return 1.17e-13 * pow(T4, (-0.942 - 0.031 * log(T4))); // cm^3 s^-1
}

double alphaB_Hbeta(double T4) {
    return 3.03e-14 * pow(T4, (-0.874 - 0.058 * log(T4))); // cm^3 s^-1
}

double alpha1_HeI(double T4) {
    return 1.54e-13 * pow(T4, -0.486); // cm^3 s^-1
}

double alphaA_HeI(double T4) {
    return alphaB_HeI(T4) + alpha1_HeI(T4);
}
//#Osterbrock Table A5.1
//alphaB_OIII   = 3.66*10**-12                                                               #cm^3 s^-1
//alphaB_OII    = 3.99*10**-13                                                               #cm^3 s^-1
//RR: Badnel2006 https://iopscience.iop.org/article/10.1086/508465
//DR: Badnell,N.R.1986,J.Phys.B,19,3827. 2006a
//RR and DR summarized in Aaron Smith COLT bitbucket
//https://bitbucket.org/aaron_smith/colt/src/fb0cd32aeadaedce637a2df46780b1a71a1d3864/src/rates.h
//#########################Oxygen#########################################
double alpha_RR_OIII(double T4) {
    double T = T4 * 10000;
    double ST0 = sqrt(T / 0.1602);
    double ST1 = sqrt(T / 4.377e6);
    double Bp = 0.7668 + 0.107 * exp(-139200. / T);
    return 2.096e-9 / (ST0 * pow((1 + ST0), (1 - Bp)) * pow((1 + ST1), (1 + Bp)));
}
double alpha_RR_OII(double T4) {
    double T = T4 * 10000;
    double ST0 = sqrt(T / 4.136);
    double ST1 = sqrt(T / 4.214e6);
    double Bp = 0.6109 + 0.4093 * exp(-87700. / T);
    return 6.622e-11 / (ST0 * pow((1 + ST0), (1 - Bp)) * pow((1 + ST1), (1 + Bp)));
}

double alpha_DR_OIII(double T4) {
    double T = T4 * 10000;
    return pow(T, -1.5) * (1.627e-7 * exp(-45.35 / T) + 1.262e-7 * exp(-284.7 / T) + 6.663e-7 * exp(-4166. / T) + 3.925e-6 * exp(-28770. / T) + 0.002406 * exp(-195300. / T) + 0.001146 * exp(-364600. / T));
}

double alpha_DR_OII(double T4) {
    double T = T4 * 10000;
    return pow(T, -1.5) * (5.629e-8 * exp(-5395. / T) + 2.55e-7 * exp(-17700. / T) + 0.0006173 * exp(-167100. / T) + 0.0001627 * exp(-268700. / T));
}

double alphaB_OIII(double T4) {
    return alpha_RR_OIII(T4) + alpha_DR_OIII(T4);
}

double alphaB_OII(double T4) {
    return alpha_RR_OII(T4) + alpha_DR_OII(T4);
}


double k0_OI_ct(double T4) {
    return 1.14e-9 * pow(T4, 0.4 + 0.018 * log(T4));
}

double k1_OI_ct(double T4) {
    return 3.44e-10 * pow(T4, 0.451 + 0.036 * log(T4));
}

double k2_OI_ct(double T4) {
    return 5.33e-10 * pow(T4, 0.384 + 0.024 * log(T4)) * exp(-97 / T4 / 10000);
}

double k0r_OI_ct(double T4) {
    return 8.0 / 5.0 * k0_OI_ct(T4) * exp(-229 / T4 / 10000);
}

// Photoionization cross section for HI
std::vector<double> sigma_HI(const std::vector<double>& nu) {
    // Define cross section vector
    std::vector<double> sigma(nu.size(), 0.0);

    // Constants for calculation
    const double E0 = 4.298e-1;     // eV
    const double sigma0 = 5.475e4;  // Mb
    const double ya = 3.288e1;
    const double P = 2.963;
    const double yw = 0.0;
    const double y0 = 0.0;
    const double y1 = 0.0;

    // Loop through frequencies
    for (size_t i = 0; i < nu.size(); ++i) {
        // Convert frequency to energy
        double E = h * nu[i] / eV2J;

        // Check energy range
        if (E < 13.6 || E > 5e4) {
            sigma[i] = 0.0;
        } else {
            // Calculate sigma for valid energy range
            double x = E / E0 - y0;
            double y = std::sqrt(x * x + y1 * y1);
            sigma[i] = sigma0 * ((x - 1) * (x - 1) + yw * yw) * std::pow(y, 0.5 * P - 5.5) * std::pow(1 + std::sqrt(y / ya), -P) * 1e-18; // cm^2
        }
    }

    return sigma;
}

// Constants for H-alpha and H-beta frequencies
const double nu_Halpha = 1.89 * eV2J / h; // Hz
const double nu_Hbeta = 2.55 * eV2J / h;  // Hz

// States and constants for HeI
const int g0_HeI = 1;
const int g1_HeI = 1;
const int g2_HeI = 3;
const int g3_HeI = 3;
const double A30_HeI = 1.26e-4; // s^-1

// Interpolation functions for collisional coefficients
std::vector<double> T4_grid = {0.6000, 0.8000, 1.0000, 1.5000, 2.0000, 2.5000};
std::vector<double> k31_grid = {1.95e-8, 2.45e-8, 2.60e-8, 3.05e-8, 2.55e-8, 2.68e-8};
std::vector<double> k32_grid = {2.34e-9, 3.64e-9, 5.92e-9, 7.83e-9, 9.23e-9, 9.81e-9};
std::vector<double> T4_grid_E30 = {3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25, 5.50, 5.75};
std::vector<double> Omega03_grid = {6.198e-2, 6.458e-2, 6.387e-2, 6.157e-2, 5.832e-2, 5.320e-2, 4.787e-2, 4.018e-2, 3.167e-2};

double k31_HeI(double T4) {
    auto it = std::upper_bound(T4_grid.begin(), T4_grid.end(), T4);
    int idx = std::distance(T4_grid.begin(), it);
    if (idx == 0) idx = 1;
    if (idx == T4_grid.size()) idx = T4_grid.size() - 1;
    double k31 = k31_grid[idx - 1] + (T4 - T4_grid[idx - 1]) / (T4_grid[idx] - T4_grid[idx - 1]) * (k31_grid[idx] - k31_grid[idx - 1]);
    return k31;
}

double k32_HeI(double T4) {
    auto it = std::upper_bound(T4_grid.begin(), T4_grid.end(), T4);
    int idx = std::distance(T4_grid.begin(), it);
    if (idx == 0) idx = 1;
    if (idx == T4_grid.size()) idx = T4_grid.size() - 1;
    double k32 = k32_grid[idx - 1] + (T4 - T4_grid[idx - 1]) / (T4_grid[idx] - T4_grid[idx - 1]) * (k32_grid[idx] - k32_grid[idx - 1]);
    return k32;
}

double k30_HeI(double T4) {
    auto it = std::upper_bound(T4_grid_E30.begin(), T4_grid_E30.end(), T4);
    int idx = std::distance(T4_grid_E30.begin(), it);
    if (idx == 0) idx = 1;
    if (idx == T4_grid_E30.size()) idx = T4_grid_E30.size() - 1;
    double Omega03 = Omega03_grid[idx - 1] + (T4 - T4_grid_E30[idx - 1]) / (T4_grid_E30[idx] - T4_grid_E30[idx - 1]) * (Omega03_grid[idx] - Omega03_grid[idx - 1]);
    return 8.629e-8 / std::sqrt(T4) * Omega03 / g3_HeI;
}

// Fraction of recombination radiation resulting in hydrogen ionization
double p(double ne, double T4) {
    double numerator = 0.75 * A30_HeI;
    double denominator = A30_HeI + ne * (k30_HeI(T4) + k31_HeI(T4) + k32_HeI(T4));
    double p_value = numerator / denominator + 0.25 * 2 / 3 + 0.75 * ne * k32_HeI(T4) / denominator;
    p_value += (0.75 * ne * k31_HeI(T4) / denominator + 0.25 * 1 / 3) * 0.56;
    return p_value;
}

// Photoionization cross section for HeI
std::vector<double> sigma_HeI(const std::vector<double>& nu) {
    // Define cross section vector
    std::vector<double> sigma(nu.size(), 0.0);

    // Constants for calculation
    const double E0 = 13.61;     // eV
    const double sigma0 = 949.2; // Mb
    const double ya = 1.469;
    const double P = 3.188;
    const double yw = 2.039;
    const double y0 = 0.4434;
    const double y1 = 2.136;

    // Loop through frequencies
    for (size_t i = 0; i < nu.size(); ++i) {
        // Convert frequency to energy
        double E = h * nu[i] / eV2J;

        // Check energy range
        if (E < 24.59 || E > 5e4) {
            sigma[i] = 0.0;
        } else {
            // Calculate sigma for valid energy range
            double x = E / E0 - y0;
            double y = std::sqrt(x * x + y1 * y1);
            sigma[i] = sigma0 * ((x - 1) * (x - 1) + yw * yw) * std::pow(y, 0.5 * P - 5.5) * std::pow(1 + std::sqrt(y / ya), -P) * 1e-18; // cm^2
        }
    }

    return sigma;
}

// Photoionization cross section for HeII
std::vector<double> sigma_HeII(const std::vector<double>& nu) {
    // Define cross section vector
    std::vector<double> sigma(nu.size(), 0.0);

    // Constants for calculation
    const double E0 = 1.72;        // eV
    const double sigma0 = 1.369e4; // Mb
    const double ya = 32.88;
    const double P = 2.963;
    const double yw = 0.0;
    const double y0 = 0.0;
    const double y1 = 0.0;

    // Loop through frequencies
    for (size_t i = 0; i < nu.size(); ++i) {
        // Convert frequency to energy
        double E = h * nu[i] / eV2J;

        // Check energy range
        if (E < 54.42 || E > 5e4) {
            sigma[i] = 0.0;
        } else {
            // Calculate sigma for valid energy range
            double x = E / E0 - y0;
            double y = std::sqrt(x * x + y1 * y1);
            sigma[i] = sigma0 * ((x - 1) * (x - 1) + yw * yw) * std::pow(y, 0.5 * P - 5.5) * std::pow(1 + std::sqrt(y / ya), -P) * 1e-18; // cm^2
        }
    }

    return sigma;
}
// Photoionization cross section for OI
std::vector<double> sigma_OI(const std::vector<double>& nu) {
    // Define cross section vector
    std::vector<double> sigma(nu.size(), 0.0);

    // Constants for calculation
    const double E0 = 1.240;        // eV
    const double sigma0 = 1.745e3;  // Mb
    const double ya = 3.784;
    const double P = 17.64;
    const double yw = 7.589e-2;
    const double y0 = 8.698;
    const double y1 = 1.271e-1;

    // Loop through frequencies
    for (size_t i = 0; i < nu.size(); ++i) {
        // Convert frequency to energy
        double E = h * nu[i] * 6.242e18;

        // Check energy range
        if (E < 13.62 || E > 538) {
            sigma[i] = 0.0;
        } else {
            // Calculate sigma for valid energy range
            double x = E / E0 - y0;
            double y = std::sqrt(x * x + y1 * y1);
            sigma[i] = sigma0 * ((x - 1) * (x - 1) + yw * yw) * std::pow(y, 0.5 * P - 5.5) * std::pow(1 + std::sqrt(y / ya), -P) * 1e-18; // cm^2
        }
    }

    return sigma;
}

// Other functions and constants for OII section
// Define state degeneracy constants for OII
const int g0_OII = 4;
const int g1_OII = 6;
const int g2_OII = 4;
const int g3_OII = 4;
const int g4_OII = 2;

// Define spontaneous decay rate constants for OII
const double A10_OII = 7.416e-06 + 3.382e-05;
const double A20_OII = 1.414e-04 + 2.209e-05;
const double A21_OII = 1.30e-07 + 1.49e-20;
const double A30_OII = 5.22e-02 + 2.43e-07;
const double A31_OII = 8.37e-03 + 9.07e-02;
const double A32_OII = 1.49e-02 + 3.85e-02;
const double A40_OII = 2.12e-02 + 3.72e-07;
const double A41_OII = 8.34e-03 + 5.19e-02;
const double A42_OII = 9.32e-03 + 7.74e-02;
const double A43_OII = 1.41e-10 + 4.24e-24;

// Level energy constants for OII
const double E10_OII = 38575 * kb; // J
const double E20_OII = 38604 * kb;
const double E30_OII = 58225 * kb;
const double E40_OII = 58228 * kb;
const double E21_OII = E20_OII - E10_OII;
const double E31_OII = E30_OII - E10_OII;
const double E32_OII = E30_OII - E20_OII;
const double E41_OII = E40_OII - E10_OII;
const double E42_OII = E40_OII - E20_OII;
const double E43_OII = E40_OII - E30_OII;

// Level energy frequency constants for OII
const double nu10_OII = E10_OII / h; // Hz
const double nu20_OII = E20_OII / h;
const double nu21_OII = E21_OII / h;
const double nu30_OII = E30_OII / h;
const double nu31_OII = E31_OII / h;
const double nu32_OII = E32_OII / h;
const double nu40_OII = E40_OII / h;
const double nu41_OII = E41_OII / h;
const double nu42_OII = E42_OII / h;
const double nu43_OII = E43_OII / h;

// Collisional (de-)excitation coefficients for OII
double Omega10_OII(double T4) {
    return 0.803 * pow(T4, 0.023 - 0.008 * log(T4));
}

double k10_OII(double T4) {
    return 8.629e-8 / sqrt(T4) * Omega10_OII(T4) / g1_OII; // cm^3 s^-1
}

double k01_OII(double T4, int g1_OII, int g0_OII, double E10_OII) {
    return g1_OII / g0_OII * k10_OII(T4) * exp(-E10_OII / (kb * T4 * 1e4));
}
double Omega20_OII(double T4) {
    return 0.550 * pow(T4, 0.054 - 0.004 * log(T4));
}

double k20_OII(double T4, int g2_OII) {
    return 8.629e-8 / sqrt(T4) * Omega20_OII(T4) / g2_OII; // cm^3 s^-1
}

double k02_OII(double T4, int g2_OII, int g0_OII, double E20_OII) {
    return g2_OII / g0_OII * k20_OII(T4, g2_OII) * exp(-E20_OII / (kb * T4 * 1e4));
}

// Define Omega21_OII, k21_OII, k12_OII in a similar manner
double Omega21_OII(double T4) {
    return 1.434 * pow(T4, -0.176 + 0.004 * log(T4));
}

double k21_OII(double T4, int g2_OII) {
    return 8.629e-8 / sqrt(T4) * Omega21_OII(T4) / g2_OII; // cm^3 s^-1
}

double k12_OII(double T4, int g2_OII, int g1_OII, double E21_OII) {
    return g2_OII / g1_OII * k21_OII(T4, g2_OII) * exp(-E21_OII / (kb * T4 * 1e4));
}

double Omega30_OII(double T4) {
    return 0.140 * pow(T4, 0.025 - 0.006 * log(T4));
}

double k30_OII(double T4, int g3_OII) {
    return 8.629e-8 / sqrt(T4) * Omega30_OII(T4) / g3_OII; // cm^3 s^-1
}

double k03_OII(double T4, int g3_OII, int g0_OII, double E30_OII) {
    return g3_OII / g0_OII * k30_OII(T4, g3_OII) * exp(-E30_OII / (kb * T4 * 1e4));
}
// Define Omega31_OII
double Omega31_OII(double T4) {
    return 0.349 * pow(T4, 0.060 + 0.052 * log(T4));
}

// Define k31_OII
double k31_OII(double T4, int g3_OII) {
    return 8.629e-8 / sqrt(T4) * Omega31_OII(T4) / g3_OII; // cm^3 s^-1
}

// Define k13_OII
double k13_OII(double T4, int g3_OII, int g1_OII, double E31_OII) {
    return g3_OII / g1_OII * k31_OII(T4, g3_OII) * exp(-E31_OII / (kb * T4 * 1e4));
}
// Define Omega32_OII
double Omega32_OII(double T4) {
    return 0.326 * pow(T4, 0.063 + 0.052 * log(T4));
}

// Define k32_OII
double k32_OII(double T4, int g3_OII, int g2_OII, double E32_OII) {
    return 8.629e-8 / sqrt(T4) * Omega32_OII(T4) / g3_OII; // cm^3 s^-1
}

// Define k23_OII
double k23_OII(double T4, int g3_OII, int g2_OII, double E32_OII) {
    return g3_OII / g2_OII * k32_OII(T4, g3_OII, g2_OII, E32_OII) * exp(-E32_OII / (kb * T4 * 1e4));
}

// Define Omega40_OII
double Omega40_OII(double T4) {
    return 0.283 * pow(T4, 0.023 - 0.004 * log(T4));
}

// Define k40_OII
double k40_OII(double T4, int g4_OII) {
    return 8.629e-8 / sqrt(T4) * Omega40_OII(T4) / g4_OII; // cm^3 s^-1
}

// Define k04_OII
double k04_OII(double T4, int g4_OII, int g0_OII, double E40_OII) {
    return g4_OII / g0_OII * k40_OII(T4, g4_OII) * exp(-E40_OII / (kb * T4 * 1e4));
}
// Define Omega41_OII
double Omega41_OII(double T4) {
    return 0.832 * pow(T4, 0.076 + 0.055 * log(T4));
}

// Define k41_OII
double k41_OII(double T4, int g4_OII, int g1_OII, double E41_OII) {
    return 8.629e-8 / sqrt(T4) * Omega41_OII(T4) / g4_OII; // cm^3 s^-1
}

// Define k14_OII
double k14_OII(double T4, int g4_OII, int g1_OII, double E41_OII) {
    return g4_OII / g1_OII * k41_OII(T4, g4_OII, g1_OII, E41_OII) * exp(-E41_OII / (kb * T4 * 1e4));
}

// Define Omega42_OII
double Omega42_OII(double T4) {
    return 0.485 * pow(T4, 0.059 + 0.052 * log(T4));
}

// Define k42_OII
double k42_OII(double T4, int g4_OII, int g2_OII, double E42_OII) {
    return 8.629e-8 / sqrt(T4) * Omega42_OII(T4) / g4_OII; // cm^3 s^-1
}

// Define k24_OII
double k24_OII(double T4, int g4_OII, int g2_OII, double E42_OII) {
    return g4_OII / g2_OII * k42_OII(T4, g4_OII, g2_OII, E42_OII) * exp(-E42_OII / (kb * T4 * 1e4));
}

// Define Omega43_OII
double Omega43_OII(double T4) {
    return 0.322 * pow(T4, 0.019 + 0.037 * log(T4));
}

// Define k43_OII
double k43_OII(double T4, int g4_OII, int g3_OII, double E43_OII) {
    return 8.629e-8 / sqrt(T4) * Omega43_OII(T4) / g4_OII; // cm^3 s^-1
}

// Define k34_OII
double k34_OII(double T4, int g4_OII, int g3_OII, double E43_OII) {
    return g4_OII / g3_OII * k43_OII(T4, g4_OII, g3_OII, E43_OII) * exp(-E43_OII / (kb * T4 * 1e4));
}
// Define R01_OII
double R01_OII(double ne, double T4) {
    return ne * k01_OII(T4, g1_OII, g0_OII, E10_OII);
}

// Define R02_OII
double R02_OII(double ne, double T4) {
    return ne * k02_OII(T4, g2_OII, g0_OII, E20_OII);
}

// Define R03_OII
double R03_OII(double ne, double T4) {
    return ne * k03_OII(T4, g3_OII, g0_OII, E30_OII);
}

// Define R04_OII
double R04_OII(double ne, double T4) {
    return ne * k04_OII(T4, g4_OII, g0_OII, E40_OII);
}

// Define R10_OII
double R10_OII(double ne, double T4) {
    return ne * k10_OII(T4) + A10_OII;
}

// Define R12_OII
double R12_OII(double ne, double T4) {
    return ne * k12_OII(T4, g2_OII, g1_OII, E21_OII);
}

// Define R13_OII
double R13_OII(double ne, double T4) {
    return ne * k13_OII(T4, g3_OII, g1_OII, E31_OII);
}

// Define R14_OII
double R14_OII(double ne, double T4) {
    return ne * k14_OII(T4, g4_OII, g1_OII, E41_OII);
}

// Define R20_OII
double R20_OII(double ne, double T4) {
    return ne * k20_OII(T4, g2_OII) + A20_OII;
}

// Define R21_OII
double R21_OII(double ne, double T4) {
    return ne * k21_OII(T4, g2_OII) + A21_OII;
}

// Define R23_OII
double R23_OII(double ne, double T4) {
    return ne * k23_OII(T4, g3_OII, g2_OII, E32_OII);
}

// Define R24_OII
double R24_OII(double ne, double T4) {
    return ne * k24_OII(T4, g4_OII, g2_OII, E42_OII);
}

// Define R30_OII
double R30_OII(double ne, double T4) {
    return ne * k30_OII(T4, g3_OII) + A30_OII;
}

// Define R31_OII
double R31_OII(double ne, double T4) {
    return ne * k31_OII(T4, g3_OII) + A31_OII;
}

// Define R32_OII
double R32_OII(double ne, double T4) {
    return ne * k32_OII(T4, g3_OII, g2_OII, E32_OII) + A32_OII;
}

// Define R34_OII
double R34_OII(double ne, double T4) {
    return ne * k34_OII(T4, g4_OII, g3_OII, E43_OII);
}

// Define R40_OII
double R40_OII(double ne, double T4) {
    return ne * k40_OII(T4, g4_OII) + A40_OII;
}

// Define R41_OII
double R41_OII(double ne, double T4) {
    return ne * k41_OII(T4, g4_OII, g1_OII, E41_OII) + A41_OII;
}

// Define R42_OII
double R42_OII(double ne, double T4) {
    return ne * k42_OII(T4, g4_OII, g2_OII, E42_OII) + A42_OII;
}

// Define R43_OII
double R43_OII(double ne, double T4) {
    return ne * k43_OII(T4, g4_OII, g3_OII, E43_OII) + A43_OII;
}

//Photoionization cross section
// Define sigma_OII function
std::vector<double> sigma_OII(const std::vector<double>& nu) {
    // Define cross section vector
    std::vector<double> sigma(nu.size(), 0.0);

    // Define energy range
    const double E_min = 35.12 * h * 6.242e18; // eV
    const double E_max = 558.1 * h * 6.242e18; // eV

    // Iterate over frequency array
    for (size_t i = 0; i < nu.size(); ++i) {
        double E = h * nu[i] * 6.242e18; // eV
        // Check if energy is within valid range
        if (E >= E_min && E <= E_max) {
            // Calculate cross section
            double E0 = 1.386; // eV
            double sigma0 = 5.967 * 10; // Mb
            double ya = 3.175 * 10;
            double P = 8.943;
            double yw = 1.934e-2;
            double y0 = 2.131 * 10;
            double y1 = 1.503e-2;

            double x = E / E0 - y0;
            double y_val = std::sqrt(x * x + y1 * y1);

            // Calculate cross section using given formula
            sigma[i] = sigma0 * ((x - 1) * (x - 1) + yw * yw) * std::pow(y_val, 0.5 * P - 5.5) * std::pow(1 + std::sqrt(y_val / ya), -P) * 1e-18;
        }
    }

    return sigma;
}

// OIII data
// State degeneracy
const int g0_OIII = 1;
const int g1_OIII = 3;
const int g2_OIII = 5;
const int g3_OIII = 5;
const int g4_OIII = 1;

// Spontaneous decay rate (s^-1)
const double A10_OIII = 2.6e-5;
const double A20_OIII = 3.5e-11;
const double A21_OIII = 9.8e-5;
const double A30_OIII = 1.9e-6;
const double A31_OIII = 0.0071;
const double A32_OIII = 0.021;
const double A40_OIII = 0;
const double A41_OIII = 0.23;
const double A42_OIII = 7.1e-4;
const double A43_OIII = 1.6;

// Level energy and frequency (in J and Hz respectively)
const double E10_OIII = 163 * kb;
const double E20_OIII = 441 * kb;
const double E30_OIII = 29169 * kb;
const double E40_OIII = 61207 * kb;
const double E21_OIII = E20_OIII - E10_OIII;
const double E31_OIII = E30_OIII - E10_OIII;
const double E32_OIII = E30_OIII - E20_OIII;
const double E41_OIII = E40_OIII - E10_OIII;
const double E42_OIII = E40_OIII - E20_OIII;
const double E43_OIII = E40_OIII - E30_OIII;
const double nu10_OIII = E10_OIII / h;
const double nu20_OIII = E20_OIII / h;
const double nu30_OIII = E30_OIII / h;
const double nu40_OIII = E40_OIII / h;
const double nu21_OIII = E21_OIII / h;
const double nu31_OIII = E31_OIII / h;
const double nu32_OIII = E32_OIII / h;
const double nu41_OIII = E41_OIII / h;
const double nu42_OIII = E42_OIII / h;
const double nu43_OIII = E43_OIII / h;

// OIII collisional (de-)excitation coefficients
// Omega and k functions
auto Omega10_OIII = [](double T4) { return 0.522 * pow(T4, (0.033 - 0.009 * log(T4))); };
auto k10_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega10_OIII(T4) / g1_OIII; };
auto k01_OIII = [](double T4) { return g1_OIII / g0_OIII * k10_OIII(T4) * exp(-E10_OIII / (kb * T4) / 10000); };

auto Omega20_OIII = [](double T4) { return 0.257 * pow(T4, (0.081 + 0.017 * log(T4))); };
auto k20_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega20_OIII(T4) / g2_OIII; };
auto k02_OIII = [](double T4) { return g2_OIII / g0_OIII * k20_OIII(T4) * exp(-E20_OIII / (kb * T4) / 10000); };

auto Omega21_OIII = [](double T4) { return 1.23 * pow(T4, (0.053 + 0.007 * log(T4))); };
auto k21_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega21_OIII(T4) / g2_OIII; };
auto k12_OIII = [](double T4) { return g2_OIII / g1_OIII * k21_OIII(T4) * exp(-E21_OIII / (kb * T4) / 10000); };

auto Omega30_OIII = [](double T4) { return 0.243 * pow(T4, (0.12 + 0.031 * log(T4))); };
auto k30_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega30_OIII(T4) / g3_OIII; };
auto k03_OIII = [](double T4) { return g3_OIII / g0_OIII * k30_OIII(T4) * exp(-E30_OIII / (kb * T4) / 10000); };

auto Omega31_OIII = [](double T4) { return 0.243 * pow(T4, (0.12 + 0.031 * log(T4))) * 3; };
auto k31_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega31_OIII(T4) / g3_OIII; };
auto k13_OIII = [](double T4) { return g3_OIII / g1_OIII * k31_OIII(T4) * exp(-E31_OIII / (kb * T4) / 10000); };

auto Omega32_OIII = [](double T4) { return 0.243 * pow(T4, (0.12 + 0.031 * log(T4))) * 5; };
auto k32_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega32_OIII(T4) / g3_OIII; };
auto k23_OIII = [](double T4) { return g3_OIII / g2_OIII * k32_OIII(T4) * exp(-E32_OIII / (kb * T4) / 10000); };

auto Omega40_OIII = [](double T4) { return 0.0321 * pow(T4, (0.118 + 0.057 * log(T4))); };
auto k40_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega40_OIII(T4) / g4_OIII; };
auto k04_OIII = [](double T4) { return g4_OIII / g0_OIII * k40_OIII(T4) * exp(-E40_OIII / (kb * T4) / 10000); };

auto Omega41_OIII = [](double T4) { return 0.0321 * pow(T4, (0.118 + 0.057 * log(T4))) * 3; };
auto k41_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega41_OIII(T4) / g4_OIII; };
auto k14_OIII = [](double T4) { return g4_OIII / g1_OIII * k41_OIII(T4) * exp(-E41_OIII / (kb * T4) / 10000); };

auto Omega42_OIII = [](double T4) { return 0.0321 * pow(T4, (0.118 + 0.057 * log(T4))) * 5; };
auto k42_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega42_OIII(T4) / g4_OIII; };
auto k24_OIII = [](double T4) { return g4_OIII / g2_OIII * k42_OIII(T4) * exp(-E42_OIII / (kb * T4) / 10000); };

auto Omega43_OIII = [](double T4) { return 0.523 * pow(T4, (0.210 - 0.099 * log(T4))); };
auto k43_OIII = [](double T4) { return 8.629e-8 / sqrt(T4) * Omega43_OIII(T4) / g4_OIII; };
auto k34_OIII = [](double T4) { return g4_OIII / g3_OIII * k43_OIII(T4) * exp(-E43_OIII / (kb * T4) / 10000); };


// Five level rates for OIII
auto R01_OIII = [](double ne, double T4) { return ne * k01_OIII(T4); };
auto R02_OIII = [](double ne, double T4) { return ne * k02_OIII(T4); };
auto R03_OIII = [](double ne, double T4) { return ne * k03_OIII(T4); };
auto R04_OIII = [](double ne, double T4) { return ne * k04_OIII(T4); };

auto R10_OIII = [](double ne, double T4) { return ne * k10_OIII(T4) + A10_OIII; };
auto R12_OIII = [](double ne, double T4) { return ne * k12_OIII(T4); };
auto R13_OIII = [](double ne, double T4) { return ne * k13_OIII(T4); };
auto R14_OIII = [](double ne, double T4) { return ne * k14_OIII(T4); };
auto R20_OIII = [](double ne, double T4) { return ne * k20_OIII(T4) + A20_OIII; };
auto R21_OIII = [](double ne, double T4) { return ne * k21_OIII(T4) + A21_OIII; };
auto R23_OIII = [](double ne, double T4) { return ne * k23_OIII(T4); };
auto R24_OIII = [](double ne, double T4) { return ne * k24_OIII(T4); };
auto R30_OIII = [](double ne, double T4) { return ne * k30_OIII(T4) + A30_OIII; };
auto R31_OIII = [](double ne, double T4) { return ne * k31_OIII(T4) + A31_OIII; };
auto R32_OIII = [](double ne, double T4) { return ne * k32_OIII(T4) + A32_OIII; };
auto R34_OIII = [](double ne, double T4) { return ne * k34_OIII(T4); };
auto R40_OIII = [](double ne, double T4) { return ne * k40_OIII(T4) + A40_OIII; };
auto R41_OIII = [](double ne, double T4) { return ne * k41_OIII(T4) + A41_OIII; };
auto R42_OIII = [](double ne, double T4) { return ne * k42_OIII(T4) + A42_OIII; };
auto R43_OIII = [](double ne, double T4) { return ne * k43_OIII(T4) + A43_OIII; };



// NII parameters
// State degeneracy
int g0_NII = 1;
int g1_NII = 3;
int g2_NII = 5;
int g3_NII = 5;
int g4_NII = 1;

// Spontaneous decay rates (s^-1)
double A10_NII = 2.08e-6;
double A20_NII = 1.12e-12;
double A21_NII = 7.46e-6;
double A30_NII = 5.25e-7;
double A31_NII = 9.22e-7 + 9.84e-4;
double A32_NII = 8.65e-6 + 2.91e-3;
double A40_NII = 0;
double A41_NII = 3.18e-2;
double A42_NII = 1.55e-4;
double A43_NII = 1.14;

// Level energy and frequency (J and s^-1)
double E10_NII = 70 * kb;
double E20_NII = 188 * kb;
double E30_NII = 22037 * kb;
double E40_NII = 47033 * kb;
double E21_NII = E20_NII - E10_NII;
double E31_NII = E30_NII - E10_NII;
double E32_NII = E30_NII - E20_NII;
double E41_NII = E40_NII - E10_NII;
double E42_NII = E40_NII - E20_NII;
double E43_NII = E40_NII - E30_NII;
double nu10_NII = E10_NII / h;
double nu20_NII = E20_NII / h;
double nu30_NII = E30_NII / h;
double nu40_NII = E40_NII / h;
double nu21_NII = E21_NII / h;
double nu31_NII = E31_NII / h;
double nu32_NII = E32_NII / h;
double nu41_NII = E41_NII / h;
double nu42_NII = E42_NII / h;
double nu43_NII = E43_NII / h;

// Collisional (de-)excitation coefficients for NII
// Omega functions
double Omega10_NII(double T4) { return 0.431 * pow(T4, 0.099 + 0.014 * log(T4)); }
double Omega20_NII(double T4) { return 0.273 * pow(T4, 0.166 + 0.030 * log(T4)); }
double Omega21_NII(double T4) { return 1.15 * pow(T4, 0.137 + 0.024 * log(T4)); }
double Omega30_NII(double T4) { return 0.303 * pow(T4, 0.053 + 0.009 * log(T4)); }
double Omega31_NII(double T4) { return 0.909 * pow(T4, 0.053 + 0.010 * log(T4)); }
double Omega32_NII(double T4) { return 1.51 * pow(T4, 0.054 + 0.011 * log(T4)); }
double Omega40_NII(double T4) { return 0.0352 * pow(T4, 0.066 + 0.018 * log(T4)); }
double Omega41_NII(double T4) { return 0.105 * pow(T4, 0.070 + 0.021 * log(T4)); }
double Omega42_NII(double T4) { return 0.176 * pow(T4, 0.065 + 0.017 * log(T4)); }
double Omega43_NII(double T4) { return 0.806 * pow(T4, -0.175 - 0.014 * log(T4)); }

// Rate coefficients
double k10_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega10_NII(T4) / g1_NII; }
double k01_NII(double T4) { return g1_NII / g0_NII * k10_NII(T4) * exp(-E10_NII / (kb * T4) / 10000); }
double k20_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega20_NII(T4) / g2_NII; }
double k02_NII(double T4) { return g2_NII / g0_NII * k20_NII(T4) * exp(-E20_NII / (kb * T4) / 10000); }
double k21_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega21_NII(T4) / g2_NII; }
double k12_NII(double T4) { return g2_NII / g1_NII * k21_NII(T4) * exp(-E21_NII / (kb * T4) / 10000); }
double k30_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega30_NII(T4) / g3_NII; }
double k03_NII(double T4) { return g3_NII / g0_NII * k30_NII(T4) * exp(-E30_NII / (kb * T4) / 10000); }
double k31_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega31_NII(T4) / g3_NII; }
double k13_NII(double T4) { return g3_NII / g1_NII * k31_NII(T4) * exp(-E31_NII / (kb * T4) / 10000); }
double k32_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega32_NII(T4) / g3_NII; }
double k23_NII(double T4) { return g3_NII / g2_NII * k32_NII(T4) * exp(-E32_NII / (kb * T4) / 10000); }
double k40_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega40_NII(T4) / g4_NII; }
double k04_NII(double T4) { return g4_NII / g0_NII * k40_NII(T4) * exp(-E40_NII / (kb * T4) / 10000); }
double k41_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega41_NII(T4) / g4_NII; }
double k14_NII(double T4) { return g4_NII / g1_NII * k41_NII(T4) * exp(-E41_NII / (kb * T4) / 10000); }
double k42_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega42_NII(T4) / g4_NII; }
double k24_NII(double T4) { return g4_NII / g2_NII * k42_NII(T4) * exp(-E42_NII / (kb * T4) / 10000); }
double k43_NII(double T4) { return 8.629e-8 / sqrt(T4) * Omega43_NII(T4) / g4_NII; }
double k34_NII(double T4) { return g4_NII / g3_NII * k43_NII(T4) * exp(-E43_NII / (kb * T4) / 10000); }

// Five-level rates for NII
double R01_NII(double ne, double T4) { return ne * k01_NII(T4); }
double R02_NII(double ne, double T4) { return ne * k02_NII(T4); }
double R03_NII(double ne, double T4) { return ne * k03_NII(T4); }
double R04_NII(double ne, double T4) { return ne * k04_NII(T4); }
double R10_NII(double ne, double T4) { return ne * k10_NII(T4) + A10_NII; }
double R12_NII(double ne, double T4) { return ne * k12_NII(T4); }
double R13_NII(double ne, double T4) { return ne * k13_NII(T4); }
double R14_NII(double ne, double T4) { return ne * k14_NII(T4); }
double R20_NII(double ne, double T4) { return ne * k20_NII(T4) + A20_NII; }
double R21_NII(double ne, double T4) { return ne * k21_NII(T4) + A21_NII; }
double R23_NII(double ne, double T4) { return ne * k23_NII(T4); }
double R24_NII(double ne, double T4) { return ne * k24_NII(T4); }
double R30_NII(double ne, double T4) { return ne * k30_NII(T4) + A30_NII; }
double R31_NII(double ne, double T4) { return ne * k31_NII(T4) + A31_NII; }
double R32_NII(double ne, double T4) { return ne * k32_NII(T4) + A32_NII; }
double R34_NII(double ne, double T4) { return ne * k34_NII(T4); }
double R40_NII(double ne, double T4) { return ne * k40_NII(T4) + A40_NII; }
double R41_NII(double ne, double T4) { return ne * k41_NII(T4) + A41_NII; }
double R42_NII(double ne, double T4) { return ne * k42_NII(T4) + A42_NII; }
double R43_NII(double ne, double T4) { return ne * k43_NII(T4) + A43_NII; }

// Cross section for NI
std::vector<double> sigma_NI(const std::vector<double>& nu) {
    std::vector<double> sigma(nu.size(), 0.0);

    const double E0 = 4.034;  // eV
    const double sigma0 = 8.235 * 100;  // Mb
    const double ya = 8.033 * 10;
    const double P = 3.928;
    const double yw = 9.097 * pow(10, -2);
    const double y0 = 8.598 * pow(10, -1);
    const double y1 = 2.325;

    for (size_t i = 0; i < nu.size(); ++i) {
        double E = h * nu[i] * 6.242 * pow(10, 18);
        if (E >= 14.53 && E <= 404.8) {
            double x = E / E0 - y0;
            double y = sqrt(pow(x, 2) + pow(y1, 2));
            sigma[i] = sigma0 * (pow(x - 1, 2) + pow(yw, 2)) * pow(y, 0.5 * P - 5.5) * pow(1 + sqrt(y / ya), -P) * pow(10, -18);  // cm^2
        }
    }

    return sigma;
}

// Cross section for NII
std::vector<double> sigma_NII(const std::vector<double>& nu) {
    std::vector<double> sigma(nu.size(), 0.0);

    const double E0_NII = 6.128 * pow(10, -2);  // eV
    const double sigma0_NII = 1.944;            // Mb
    const double ya_NII = 8.163 * pow(10, 2);
    const double P_NII = 8.773;
    const double yw_NII = 1.043 * pow(10, 1);
    const double y0_NII = 4.280 * pow(10, 2);
    const double y1_NII = 2.030 * pow(10, 1);

    for (size_t i = 0; i < nu.size(); ++i) {
        double E_NII = h * nu[i] * 6.242 * pow(10, 18);
        if (E_NII >= 29.6 && E_NII <= 423.6) {
            double x_NII = E_NII / E0_NII - y0_NII;
            double y_NII = sqrt(pow(x_NII, 2) + pow(y1_NII, 2));
            sigma[i] = sigma0_NII * (pow(x_NII - 1, 2) + pow(yw_NII, 2)) * pow(y_NII, 0.5 * P_NII - 5.5) * pow(1 + sqrt(y_NII / ya_NII), -P_NII) * pow(10, -18);  // cm^2
        }
    }

    return sigma;
}

// Cross section for NIII
std::vector<double> sigma_NIII(const std::vector<double>& nu) {
    std::vector<double> sigma(nu.size(), 0.0);

    const double E0_NIII = 0.2420;  // eV
    const double sigma0_NIII = 0.9375;  // Mb
    const double ya_NIII = 278.8;
    const double P_NIII = 9.156;
    const double yw_NIII = 1.850;
    const double y0_NIII = 187.7;
    const double y1_NIII = 3.999;

    for (size_t i = 0; i < nu.size(); ++i) {
        double E_NIII = h * nu[i] * 6.242 * pow(10, 18);
        if (E_NIII >= 47.45 && E_NIII <= 447.3) {
            double x_NIII = E_NIII / E0_NIII - y0_NIII;
            double y_NIII = sqrt(pow(x_NIII, 2) + pow(y1_NIII, 2));
            sigma[i] = sigma0_NIII * (pow(x_NIII - 1, 2) + pow(yw_NIII, 2)) * pow(y_NIII, 0.5 * P_NIII - 5.5) * pow(1 + sqrt(y_NIII / ya_NIII), -P_NIII) * pow(10, -18);  // cm^2
        }
    }

    return sigma;
}