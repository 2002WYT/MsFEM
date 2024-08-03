#include "head.h"
#include <cmath>
const double eps = 1.0/32;
#define pi M_PI
double get_kappa(double x, double y, const int example) {
    switch (example)
    {
        case 1: {return 3 + x*y; break;}
        case 2: {return (2.0+1.8*sin(2*M_PI*x/eps))/(2.0+1.8*cos(2*M_PI*y/eps)) + 
                (2.0+sin(2*M_PI*y/eps))/(2.0+1.8*sin(2*M_PI*x/eps)); break;}
        case 3: {return 1.8 + sin(2*M_PI*(x+y)/eps); break;}
        case 4: {return (2.0+sin(x/eps))/(2+1.8*cos(y/eps)); break;}
        case 5: {return 1.0/(1.5+sin(2*pi*x/eps))+1.0/(1.5+cos(2*pi*y/eps)); break;}
        case 6: 
        {return (1.3+sin(2*pi*x/eps))*(1.3+sin(2*pi*y/eps))*(1.3+cos(2*pi*x/eps+2*pi/5))*
            (1.3+cos(2*pi*y/eps+2*pi/5)); break;}
        default: return -1;
    }
    
}