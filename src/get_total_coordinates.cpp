#include "head.h"
void get_total_coordinates_x(const int nx, const int ny, const int gs, const double x0,
    const double hx, const double hy, 
    const double* pol_x, const double* pol_y, double* coorx)
{
    for(int i = 0; i < gs*gs*ny; i++) {
        for(int j = 0; j < nx; j++) {
            coorx[j*gs*gs*ny+i] = pol_x[i%(gs*gs)/gs] + hx/nx*j + x0;
        }
    }
}
void get_total_coordinates_y(const int nx, const int ny, const int gs, const double y0, 
    const double hx, const double hy, 
    const double* pol_x, const double* pol_y, double* coory)
{
    for(int i = 0; i < gs*gs*ny; i++) {
        for(int j = 0; j < nx; j++) {
            coory[j*gs*gs*ny+i] = pol_y[i%gs] + (int)(i/(gs*gs))*hy/ny + y0;
        }
    }
}