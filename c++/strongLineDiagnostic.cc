#include <iostream>
#include <vector>
#include <cmath>
//#include "const_ion.h"
#include "HIILines.cc"

using namespace std;

vector<double> MeasurelogZ(double R_Obs, double O3O2_Obs, double Rmode, string tag) {
    vector<double> result;
    vector<double> logZQ(1000);
    for (int i = 0; i < 1000; i++) {
        logZQ[i] = -2 + i * 0.002;
    }
    vector<double> ZQ(1000);
    for (int i = 0; i < 1000; i++) {
        ZQ[i] = pow(10, logZQ[i]);
    }
    vector<double> T4OIII(1000);
    vector<double> T4OII(1000);
    vector<int> idx_sol;
    while (idx_sol.size() == 0) {
            double fO3O2 = 1 / (0.75 * k03_OIII(T4OIII) * nu32_OIII / (k01_OII(T4OII) * nu10_OII + k02_OII(T4OII) * nu20_OII));
            double FOIII = fO3O2 * O3O2_Obs / (1 + fO3O2 * O3O2_Obs);
            double FOII = 1 - FOIII;
            double T4HII = T4OIII * FOIII + T4OII * FOII;
            double fR2 = pow(10, -3.31) / nu_Hbeta / alphaB_Hbeta(T4HII) * (k01_OII(T4OII) * nu10_OII + k02_OII(T4OII) * nu20_OII);
            double R2 = fR2 * FOII * ZQ;
            double fR3 = 0.75 * pow(10, -3.31) * k03_OIII(T4OIII) * nu32_OIII / nu_Hbeta / alphaB_Hbeta(T4HII);
            double R3 = fR3 * FOIII * ZQ;
            double fR3p = 0.25 * pow(10, -3.31) * k03_OIII(T4OIII) * nu31_OIII / nu_Hbeta / alphaB_Hbeta(T4HII);
            double R3p = fR3p * FOIII * ZQ;
            double R23 = R2 + R3 + R3p;
            if (Rmode == “R2”) {
                R = R2;
            } else if (Rmode == “R3”) {
                R = R3;
            } else if (Rmode == “R23”) {
                R = R23;
            }
            if (*max_element(R.begin(), R.end()) >= R_Obs) {
                vector<double> diff;
                for (int i = 0; i < R.size() - 1; i++) {
                    diff.push_back(R[i + 1] - R[i]);
                }
                vector<double> diff_prod;
                for (int i = 0; i < diff.size() - 1; i++) {
                    diff_prod.push_back(diff[i] * diff[i + 1]);
                }
                for (int i = 0; i < diff_prod.size(); i++) {
                    if (diff_prod[i] <= 0) {
                        idx_sol.push_back(i);
                        break;
                    }
                }
            }
            if (*max_element(R.begin(), R.end()) < R_Obs) {
                T4OIII += 0.1;
                T4OII += 0.1;
            }
        }
        double w1 = (R[idx_sol[0] + 1] - R_Obs) / (R[idx_sol[0] + 1] - R[idx_sol[0]);
        double w2 = 1 - w1;
        
        vector<double> result;
        result.push_back(w1 * logZQ[idx_sol[0]] + w2 * logZQ[idx_sol[0] + 1]);
        return result;
    }

    vector<double> MeasurelogZ_withN2O2(double R_Obs, double O3O2_Obs, double N2O2, double Rmode, string tag) {
    vector<double> result;
    vector<double> logZQ(1000);
    for (int i = 0; i < 1000; i++) {
        logZQ[i] = -2 + i * 0.002;
    }
    vector<double> ZQ(1000);
    for (int i = 0; i < 1000; i++) {
        ZQ[i] = pow(10, logZQ[i]);
    }
    vector<double> T4OIII(1000);
    vector<double> T4OII(1000);
    // Implement N2O2 model calculation
    vector<double> N2O2_model;
    for (int i = 0; i < 1000; i++) {
        double T4OII_sol = w1T4OII[i] + w2T4OII[i + 1];
        double T4OIII_sol = w1T4OIII[i] + w2T4OIII[i + 1];
        double temp_N2O2_model = 0.0738 * k03_NII(T4OII_sol) / (k01_OII(T4OII_sol) + k02_OII(T4OII_sol));
        N2O2_model.push_back(temp_N2O2_model);
    }
    // Find the index where the model N2O2 is closest to the observed N2O2
    double min_diff = abs(N2O2_model[0] - N2O2);
    int idx_min_diff = 0;
    for (int i = 1; i < 1000; i++) {
        double temp_diff = abs(N2O2_model[i] - N2O2);
        if (temp_diff < min_diff) {
            min_diff = temp_diff;
            idx_min_diff = i;
        }
    }
    // Interpolate logZQ based on the closest N2O2 model
    double w1_N2O2 = (N2O2_model[idx_min_diff] - N2O2) / (N2O2_model[idx_min_diff] - N2O2_model[idx_min_diff - 1]);
    double w2_N2O2 = 1 - w1_N2O2;
    result.push_back(w1_N2O2 * logZQ[idx_min_diff] + w2_N2O2 * logZQ[idx_min_diff - 1]);
    return result;
}


//###########################################

int main() {
    vector<double> logZQ = MeasurelogZ(10.0, 0.5, “R2”, “local”);
    vector<double> logZQN2 = MeasurelogZ_withN2O2(10.0, 0.5, 0.4, “R3”, “EoR”);

    cout << “Result of MeasurelogZ function:” << endl;
    for (int i = 0; i < logZQ.size(); i++) {
        cout << logZQ[i] << " ";
    }
    cout << endl;

    cout << “Result of MeasurelogZ_withN2O2 function:” << endl;
    for (int i = 0; i < logZQN2.size(); i++) {
        cout << logZQN2[i] << " ";
    }
    cout << endl;
    return 0;
}