#include "head.h"
void get_real_solution(double* meshx, double* meshy, double* ue, int sizemeshgrid, 
    const int example) {
    for (int i = 0; i < sizemeshgrid; i++) {
        ue[i] = real_u(meshx[i], meshy[i], example); 
    }
}
