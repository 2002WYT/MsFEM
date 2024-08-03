#include "head.h"
void get_meshgrid(const int nx, const int ny, const int order, const double x0, const double y0, 
    const double xl, const double yl, double* meshx, double* meshy)
{
    for (int i = 0; i < nx*order+1; i++) {
        for (int j = 0; j < ny*order+1; j++) {
            meshx[i*(ny*order+1)+j] = x0 + i*(xl-x0)/(nx*order);
            meshy[i*(ny*order+1)+j] = y0 + j*(yl-y0)/(ny*order);
        }
    }
}