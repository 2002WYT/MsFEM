#include "head.h"
#include <cmath>
//#include <iostream>
double get_relative_L2_norm(double* u, double* ue, int sizeu) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < sizeu; i++) {
        sum1 += pow(u[i] - ue[i], 2);
        sum2 += pow(ue[i], 2);
    }
    //std::cout << "sum: " << sum << std::endl; 
    return sqrt(sum1/sum2);
}