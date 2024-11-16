// #include "const.h"
// #include "const_ion.h"
#include "func_declarations.h"
#include <iostream>
#include "HIILines.cc"
// #include "solve_levels.cc"

int main(int argc, char* argv[]) {
    
    // Check if temperature argument T4 is provided
    if (argc < 2) {
        auto someFunction = [](double ne, double T4) { /* ... */ };
        return 1;
    }

    // Convert command-line string to double (temperature scaled by 10,000 K)
    std::string temp_str(argv[1]);
    double T4 = std::stod(temp_str);

    // Calculate and print the recombination coefficients for different species at temperature T4
    std::cout << "alphaB_HI(T4): " << alphaB_HI(T4) << " cm^3 s^-1" << std::endl;
    std::cout << "alphaB_HeI(T4): " << alphaB_HeI(T4) << " cm^3 s^-1" << std::endl;
    std::cout << "alphaB_Halpha(T4): " << alphaB_Halpha(T4) << " cm^3 s^-1" << std::endl;
    std::cout << "alphaB_Hbeta(T4): " << alphaB_Hbeta(T4) << " cm^3 s^-1" << std::endl;
    std::cout << "alphaB_OIII(T4): " << alphaB_OIII(T4) << " cm^3 s^-1" << std::endl;
    std::cout << "alphaB_OII(T4): " << alphaB_OII(T4) << " cm^3 s^-1" << std::endl;

    ////////////////////////////////////////////
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
    // solve_levels();

}