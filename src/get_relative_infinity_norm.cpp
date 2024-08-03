#include "head.h"
#include <cmath>
//#include <iostream>
double get_relative_infinity_norm(double* u, double* ue, int sizeu) {
    double maxi = 0.0;
    double maxe = 0.0;
    for (int i = 0; i < sizeu; i++)
    {
        if (fabs(u[i] - ue[i]) > maxi) maxi = fabs(u[i] - ue[i]);
        if (fabs(ue[i]) > maxe) maxe = fabs(ue[i]);
    }
    return maxi/maxe;
}