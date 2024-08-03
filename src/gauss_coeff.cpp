#include "head.h"

void gauss_coefficient(const int gs, double* coeff) {
    double COE[5][5] = {
        {2.0},
        {1.0, 1.0},
        {0.5555555555555556, 0.8888888888888888, 0.5555555555555556},
        {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538},
        {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891}
    };
    for (int i = 0; i < gs; i++) {
        coeff[i] = COE[gs - 1][i];
    }
}