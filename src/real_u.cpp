#include "head.h"
#include <cmath>
#define pi M_PI
extern const double eps;
double real_u(double x, double y, const int example) {
    switch (example)
    {
        case 1: {return x*(1-x)*y*(1-y); break;}
        default: return -1;
        case 4: {return sin(100*pi*x)*sin(50*pi*y); break;}
    }
    
}