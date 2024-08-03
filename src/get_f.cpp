#include "head.h"
#include <cmath>
#define pi M_PI
extern const double eps;
double get_f(double x, double y, const int example) {
  switch (example)
    {
        case 1: {return - x*(x*y*(x - 1) + x*(x - 1)*(y - 1)) - y*(x*y*(y - 1) + y*(x - 1)*(y - 1)) - 
          2*x*(x*y + 3)*(x - 1) - 2*y*(x*y + 3)*(y - 1); break;}
        case 2: {return 1.0; break;}
        case 3: {return - 2*x*(sin((2*M_PI*(x + y))/eps) + 1.8)*(x - 1) - 2*y*(sin((2*M_PI*(x + y))
          /eps) + 1.8)*(y - 1) - (2*M_PI*cos((2*M_PI*(x + y))/eps)*(x*y*(x - 1) + x*(x - 1)*(y - 1)))
          /eps - (2*M_PI*cos((2*M_PI*(x + y))/eps)*(x*y*(y - 1) + y*(x - 1)*(y - 1)))/eps; break;}

        case 4: {
          return 12500*pi*pi*sin(100*pi*x)*sin(50*pi*y)*(sin(x/eps)*cos((2*pi)/5 + y/eps) 
          + 1.8) - (100*pi*cos((2*pi)/5 + y/eps)*cos(100*pi*x)*sin(50*pi*y)*cos(x/eps))
          /eps + (50*pi*sin(x/eps)*sin((2*pi)/5 + y/eps)*cos(50*pi*y)*sin(100*pi*x))/eps;
          break;
        case 5: {return 1.0; break;}
        case 6: {return 1.0; break;}

        }
        default: return -1;
    }
}