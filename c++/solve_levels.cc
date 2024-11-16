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

using std::cout, std::endl, std::setprecision, std::vector, std::array;

// ################# COMPLETE BASED ON THE ACTUAL DATA ###############

static const bool verbose = true; // Verbose output
// static double lognH; // Replace 1.0 with the actual value for lognH
// static double nH;
// static double logZ; // Replace -2.0 with the actual value for logZ
// static double Z;
// static double T4; // Replace 1.0 with the actual value for T4OII
// lognH = 0.0; // Use the global variable lognH
// nH = pow(10, lognH);
// logZ = -2.0; // Use the global variable logZ
// Z = pow(10, logZ);
// T4 = 1.0; // T/10^4 K

double R11_OIII(double ne, double T4) { return -R10_OIII(ne, T4) - R12_OIII(ne, T4) - R13_OIII(ne, T4) - R14_OIII(ne, T4); }
double R22_OIII(double ne, double T4) { return -R20_OIII(ne, T4) - R21_OIII(ne, T4) - R23_OIII(ne, T4) - R24_OIII(ne, T4); }
double R33_OIII(double ne, double T4) { return -R30_OIII(ne, T4) - R31_OIII(ne, T4) - R32_OIII(ne, T4) - R34_OIII(ne, T4); }
double R44_OIII(double ne, double T4) { return -R40_OIII(ne, T4) - R41_OIII(ne, T4) - R42_OIII(ne, T4) - R43_OIII(ne, T4); }

// Function pointers for rate coefficients
static constexpr int n_levels = 5; // Number of levels
constexpr int nm1 = n_levels - 1; // Number of excited levels
static vector<vector<double (*)(double, double)>> Rij_functions = {
  {R11_OIII, R21_OIII, R31_OIII, R41_OIII},
  {R12_OIII, R22_OIII, R32_OIII, R42_OIII},
  {R13_OIII, R23_OIII, R33_OIII, R43_OIII},
  {R14_OIII, R24_OIII, R34_OIII, R44_OIII}
};
static const array<double (*)(double, double), nm1> R0i_functions = {R01_OIII, R02_OIII, R03_OIII, R04_OIII};

// Function to solve OIII level population abundances
// Defaults: ne = 1 cm^-3, T4 = 1 K
// Solves A x = b, where A is the rate coefficient matrix, x is the level population vector, and b is the source vector
void solve_levels(double ne = 10000., double T4 = 1.) {
  gsl_matrix *A = gsl_matrix_alloc(nm1, nm1);
  gsl_vector *b = gsl_vector_alloc(nm1);
  gsl_vector *x = gsl_vector_alloc(nm1);

  // Set up the matrix A (rate coefficient matrix)
  for (int i = 0; i < nm1; ++i)
    for (int j = 0; j < nm1; ++j)
      gsl_matrix_set(A, i, j, -Rij_functions[i][j](ne, T4)); // Set matrix elements
  // gsl_matrix_set(A, 0, 0, R10_OIII(ne, T4) + R12_OIII(ne, T4) + R13_OIII(ne, T4) + R14_OIII(ne, T4));
  // gsl_matrix_set(A, 0, 1, -R21_OIII(ne, T4));
  // gsl_matrix_set(A, 0, 2, -R31_OIII(ne, T4));
  // gsl_matrix_set(A, 0, 3, -R41_OIII(ne, T4));
  // gsl_matrix_set(A, 1, 0, -R12_OIII(ne, T4));
  // gsl_matrix_set(A, 1, 1, R20_OIII(ne, T4) + R21_OIII(ne, T4) + R23_OIII(ne, T4) + R24_OIII(ne, T4));
  // gsl_matrix_set(A, 1, 2, -R32_OIII(ne, T4));
  // gsl_matrix_set(A, 1, 3, -R42_OIII(ne, T4));
  // gsl_matrix_set(A, 2, 0, -R13_OIII(ne, T4));
  // gsl_matrix_set(A, 2, 1, -R23_OIII(ne, T4));
  // gsl_matrix_set(A, 2, 2, R30_OIII(ne, T4) + R31_OIII(ne, T4) + R32_OIII(ne, T4) + R34_OIII(ne, T4));
  // gsl_matrix_set(A, 2, 3, -R43_OIII(ne, T4));
  // gsl_matrix_set(A, 3, 0, -R14_OIII(ne, T4));
  // gsl_matrix_set(A, 3, 1, -R24_OIII(ne, T4));
  // gsl_matrix_set(A, 3, 2, -R34_OIII(ne, T4));
  // gsl_matrix_set(A, 3, 3, R40_OIII(ne, T4) + R41_OIII(ne, T4) + R42_OIII(ne, T4) + R43_OIII(ne, T4));

  // Set up vector b (source vector, rates from ground state to other levels)
  for (int i = 0; i < nm1; ++i)
    gsl_vector_set(b, i, R0i_functions[i](ne, T4)); // Set vector elements
  // gsl_vector_set(b, 0, R01_OIII(ne, T4));
  // gsl_vector_set(b, 1, R02_OIII(ne, T4));
  // gsl_vector_set(b, 2, R03_OIII(ne, T4));
  // gsl_vector_set(b, 3, R04_OIII(ne, T4));

  // Solve for level populations
  gsl_permutation *p = gsl_permutation_alloc(4);
  int s;
  gsl_linalg_LU_decomp(A, p, &s); // LU decomposition
  gsl_linalg_LU_solve(A, p, b, x); // Solve for level populations

  // Calculate n0, n1, n2, n3, n4
  // Calculate sum of upper levels: Σ_{i=1} n_i
  double n_sum = 0.;
  for (int i = 0; i < nm1; ++i)
    n_sum += gsl_vector_get(x, i);
  vector<double> y(n_levels);
  y[0] = 1. / (1. + n_sum); // Ground state population (relative to ion density)
  for (int i = 0; i < nm1; ++i)
    y[i+1] = y[0] * gsl_vector_get(x, i); // Excited state populations (relative to ion density)

  // Print level populations
  if (verbose) {
    cout << setprecision(10);
    cout << "OIII level populations:" << endl;
    for (int i = 0; i < n_levels; ++i)
      cout << " n_" << i << " = " << y[i] << endl;
  }

  // Free allocated memory
  gsl_vector_free(x);
  gsl_vector_free(b);
  gsl_matrix_free(A);
}


// Main function
int main() {
  solve_levels();
  // Print the result
  // cout << "Result from assignL function:" << endl;
  // for (double val : result) {
  //     cout << val << " ";
  // }
  return 0;
}


