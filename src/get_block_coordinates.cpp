#include "head.h"
void get_block_coordinates_x(const int Nx, const int Ny, const int nx, const int ny, 
    const double x0, const int gs, 
    const double hx, const double hy, const double* pol_x, const double* pol_y, double* cocox)
    {
        for (int i = 0; i < Nx*Ny; i++)
        {
            double xi0 = x0 + (int)(i/Ny) * hx / Nx;
            for (int j = 0; j < nx*ny; j++)
            {
                double xj0 = xi0 + (int)(j/ny) * hx / Nx / nx;
                for (int k = 0; k < gs*gs; k++)
                {
                    cocox[i*nx*ny*gs*gs+j*gs*gs+k] = pol_x[k/gs] + xj0;
                }
            }
        }
    }
void get_block_coordinates_y(const int Nx, const int Ny, const int nx, const int ny,
    const double y0, const int gs, 
    const double hx, const double hy, const double* pol_x, const double* pol_y, double* cocoy)
    {
        for (int i = 0; i < Nx*Ny; i++)
        {
            double yi0 = y0 + (int)(i%Ny) * hy / Ny;
            for (int j = 0; j < nx*ny; j++)
            {
                double yj0 = yi0 + (int)(j%ny) * hy / Ny / ny;
                for (int k = 0; k < gs*gs; k++)
                {
                    cocoy[i*nx*ny*gs*gs+j*gs*gs+k] = pol_y[k%gs] + yj0;
                }
            }
        }
    }
