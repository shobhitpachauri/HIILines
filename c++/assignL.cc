#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_log.h>
#include "const.h"
#include "const_ion.h"
#include "func_declarations.h"

using namespace std;

// ################# COMPLETE BASED ON THE ACTUAL DATA ###############

// ################# COMPLETE BASED ON THE ACTUAL DATA ###############

double lognH = 1.0; // Replace 1.0 with the actual value for lognH
double logne = 1.0; // Replace 1.0 with the actual value for logne
double nH = pow(10, lognH);
double nu = 1.0; // Replace 1.0 with the actual value for nu
double logZ = 1.0; // Replace 1.0 with the actual value for logZ
double Z = pow(10, logZ);
double logQ= 1.0; // Replace 1.0 with the actual value for logQ
double T4OII = 1.0; // Replace 1.0 with the actual value for T4OII
double VOII2VHII = 1.0; // Replace 1.0 with the actual value for VOII2VHII
double VOIII2VHII = 1.0; // Replace 1.0 with the actual value for VOIII2VHII
double T4OIII = 1.0; // Replace 1.0 with the actual value for T4OIII
double T4HII = T4OII * VOII2VHII + T4OIII * VOIII2VHII;


// Function to solve OIII level population abundances
vector<double> assignL(std::vector<double> INPUT) {
    
    std::cout << std::setprecision(10);
    
    double logQ = INPUT[0], lognH = INPUT[1], logZ = INPUT[2];
    double T4OII = INPUT[3], T4OIII = INPUT[4], VOII2VHII = INPUT[5], VOIII2VHII = INPUT[6];
    double T4HII = T4OII * VOII2VHII + T4OIII * VOIII2VHII;
    double nO = pow(10, -3.31) * (pow(10, logZ)) * pow(10, lognH);
    double ne = pow(10, lognH) * (1 + 0.0737 + 0.0293 * (pow(10, logZ)));
    
    gsl_matrix *A = gsl_matrix_alloc(4, 4);
    gsl_vector *B = gsl_vector_alloc(4);
    gsl_vector *levelPop = gsl_vector_alloc(4);
   
    
    // Set up the matrix A
    gsl_matrix_set(A, 0, 0, R10_OIII(ne, T4OIII) + R12_OIII(ne, T4OIII) + R13_OIII(ne, T4OIII) + R14_OIII(ne, T4OIII));
    gsl_matrix_set(A, 0, 1, -R21_OIII(ne, T4OIII));
    gsl_matrix_set(A, 0, 2, -R31_OIII(ne, T4OIII));
    gsl_matrix_set(A, 0, 3, -R41_OIII(ne, T4OIII));
    gsl_matrix_set(A, 1, 0, -R12_OIII(ne, T4OIII));
    gsl_matrix_set(A, 1, 1, R20_OIII(ne, T4OIII) + R21_OIII(ne, T4OIII) + R23_OIII(ne, T4OIII) + R24_OIII(ne, T4OIII));
    gsl_matrix_set(A, 1, 2, -R32_OIII(ne, T4OIII));
    gsl_matrix_set(A, 1, 3, -R42_OIII(ne, T4OIII));
    gsl_matrix_set(A, 2, 0, -R13_OIII(ne, T4OIII));
    gsl_matrix_set(A, 2, 1, -R23_OIII(ne, T4OIII));
    gsl_matrix_set(A, 2, 2, R30_OIII(ne, T4OIII) + R31_OIII(ne, T4OIII) + R32_OIII(ne, T4OIII) + R34_OIII(ne, T4OIII));
    gsl_matrix_set(A, 2, 3, -R43_OIII(ne, T4OIII));
    gsl_matrix_set(A, 3, 0, -R14_OIII(ne, T4OIII));
    gsl_matrix_set(A, 3, 1, -R24_OIII(ne, T4OIII));
    gsl_matrix_set(A, 3, 2, -R34_OIII(ne, T4OIII));
    gsl_matrix_set(A, 3, 3, R40_OIII(ne, T4OIII) + R41_OIII(ne, T4OIII) + R42_OIII(ne, T4OIII) + R43_OIII(ne, T4OIII));

    // Set up vector B
    gsl_vector_set(B, 0, R01_OIII(ne, T4OIII));
    gsl_vector_set(B, 1, R02_OIII(ne, T4OIII));
    gsl_vector_set(B, 2, R03_OIII(ne, T4OIII));
    gsl_vector_set(B, 3, R04_OIII(ne, T4OIII));
    // Set up other elements of vector B
    // ...
    // Solve for level populations
 
    gsl_permutation *p = gsl_permutation_alloc(4);
    int s;
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, B, levelPop);

    // Calculate n0, n1, n2, n3, n4
    double n0 = nO / (1 + gsl_vector_get(levelPop, 0) + gsl_vector_get(levelPop, 1) + gsl_vector_get(levelPop, 2) + gsl_vector_get(levelPop, 3));
    double n1 = n0 * gsl_vector_get(levelPop, 0);
    double n2 = n0 * gsl_vector_get(levelPop, 1);
    double n3 = n0 * gsl_vector_get(levelPop, 2);
    double n4 = n0 * gsl_vector_get(levelPop, 3);
  
    // Calculate OIII line luminosities
    double logL10_OIII = log10(h) + log10(nu10_OIII * n1 * A10_OIII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOIII2VHII);
    double logL21_OIII = log10(h) + log10(nu21_OIII * n2 * A21_OIII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOIII2VHII);
    double logL31_OIII = log10(h) + log10(nu31_OIII * n3 * A31_OIII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOIII2VHII);
    double logL32_OIII = log10(h) + log10(nu32_OIII * n3 * A32_OIII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOIII2VHII);
    double logL43_OIII = log10(h) + log10(nu43_OIII * n4 * A43_OIII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOIII2VHII);

    
    // Solve linear equations A x = B
    //gsl_permutation *p = gsl_permutation_alloc(4);
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, B, levelPop);
   
    // Solve OII level population abundances
    gsl_matrix *A_OII = gsl_matrix_alloc(4, 4);
    gsl_vector *B_OII = gsl_vector_alloc(4);

    gsl_matrix_set(A_OII, 0, 0, R10_OII(ne, T4OII) + R12_OII(ne, T4OII) + R13_OII(ne, T4OII) + R14_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 0, 1, -R21_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 0, 2, -R31_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 0, 3, -R41_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 1, 0, -R12_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 1, 1, R20_OII(ne, T4OII) + R21_OII(ne, T4OII) + R23_OII(ne, T4OII) + R24_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 1, 2, -R32_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 1, 3, -R42_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 2, 0, -R13_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 2, 1, -R23_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 2, 2, R30_OII(ne, T4OII) + R31_OII(ne, T4OII) + R32_OII(ne, T4OII) + R34_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 2, 3, -R43_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 3, 0, -R14_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 3, 1, -R24_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 3, 2, -R34_OII(ne, T4OII));
    gsl_matrix_set(A_OII, 3, 3, R40_OII(ne, T4OII) + R41_OII(ne, T4OII) + R42_OII(ne, T4OII) + R43_OII(ne, T4OII));
    
    
    gsl_vector_set(B_OII, 0, R01_OII(ne, T4OII));
    gsl_vector_set(B_OII, 1, R02_OII(ne, T4OII));
    gsl_vector_set(B_OII, 2, R03_OII(ne, T4OII));
    gsl_vector_set(B_OII, 3, R04_OII(ne, T4OII));
    // Set up other elements of vector B_OII
 
    gsl_vector *levelPop_OII = gsl_vector_alloc(4);
    gsl_linalg_LU_solve(A_OII, p, B_OII, levelPop_OII);

    double n0_OII = nO / (1 + gsl_vector_get(levelPop_OII, 0) + gsl_vector_get(levelPop_OII, 1) + gsl_vector_get(levelPop_OII, 2) + gsl_vector_get(levelPop_OII, 3));
    double n1_OII = n0_OII * gsl_vector_get(levelPop_OII, 0);
    double n2_OII = n0_OII * gsl_vector_get(levelPop_OII, 1);
    double n3_OII = n0_OII * gsl_vector_get(levelPop_OII, 2);
    double n4_OII = n0_OII * gsl_vector_get(levelPop_OII, 3);

    

    double logL10_OII = log10(h) + log10(nu10_OII * n1_OII * A10_OII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOII2VHII);
    double logL20_OII = log10(h) + log10(nu20_OII * n2_OII * A20_OII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOII2VHII);
    double logL30_OII = log10(h) + log10(nu30_OII * n3_OII * A30_OII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOII2VHII);
    double logL31_OII = log10(h) + log10(nu31_OII * n3_OII * A31_OII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOII2VHII);
    double logL32_OII = log10(h) + log10(nu32_OII * n3_OII * A32_OII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOII2VHII);
    double logL40_OII = log10(h) + log10(nu40_OII * n4_OII * A40_OII / alphaB_HI(T4HII) * Jps2Lsun) + logQ - lognH - logne + log10(VOII2VHII);

    double logLHalpha = log10(h) + log10(nu_Halpha * alphaB_Halpha(T4HII) / alphaB_HI(T4HII) * Jps2Lsun) + logQ;
    double logLHbeta = log10(h) + log10(nu_Hbeta * alphaB_Hbeta(T4HII) / alphaB_HI(T4HII) * Jps2Lsun) + logQ;
   
    // Free allocated memory
    gsl_matrix_free(A);
    gsl_vector_free(B);
    gsl_vector_free(levelPop);

    // Return the calculated line luminosities
    return {logL10_OII, logL20_OII, logL30_OII, logL31_OII, logL32_OII, logL40_OII, logL10_OIII, logL21_OIII, logL31_OIII, logL32_OIII, logL43_OIII, logLHalpha, logLHbeta};
}



// Main function

int main() {
    
    // Define the input parameters
    logQ = 1.0; // Use the global variable logQ
    lognH = 1.0; // Use the global variable lognH
    logZ = 1.0; // Use the global variable logZ
    T4OII = 1.0; // Use the global variable T4OII
    T4OIII = 1.0; // Use the global variable T4OIII
    VOII2VHII = 1.0; // Use the global variable VOII2VHII
    VOIII2VHII = 1.0; // Use the global variable VOIII2VHII
    
    // Call assignL
    std::vector<double> input_values = {logQ, lognH, logZ, T4OII, T4OIII, VOII2VHII, VOIII2VHII};
    std::vector<double> result = assignL(input_values);
    // Print the result
    std::cout << "Result from assignL function:" << std::endl;
    for (double val : result) {
        std::cout << val << " ";
    }
    return 0;
}


