
#include "head.h"
#include <iostream>
using namespace std;
void get_poly_term(const int order, const double* inpol, 
    const int nx, const int ny, const int flag, double* poly_term)
{
    if (flag < 0 || flag >= (order+1)*(order+1)) // 从0开始计数
    {
        cout << "flag is out of range" << endl;
        return;
    }
    const int gs = order + 1;
    double* pol_x1 = new double[order+1];
    double* pol_y1 = new double[order+1];
    for (int i = 0; i < order+1; i++) {
        pol_x1[i] = inpol[i];
        pol_y1[i] = inpol[i];
    }
    add_constant_to_vector(pol_x1, gs, 1.0);
    add_constant_to_vector(pol_y1, gs, 1.0);
    cblas_dscal(gs, 1.0/2*1.0/nx, (double *)pol_x1, 1);
    cblas_dscal(gs, 1.0/2*1.0/ny, (double *)pol_y1, 1);
    double* x_term = new double[gs*gs*nx*ny];
    double* y_term = new double[gs*gs*nx*ny];
    get_total_coordinates_x(nx, ny, gs, 0.0, 1.0, 1.0, pol_x1, 
        pol_y1, x_term);
    get_total_coordinates_y(nx, ny, gs, 0.0, 1.0, 1.0, pol_x1, 
        pol_y1, y_term);
    double* tx = new double[nx*order+1];
    for (int i = 0; i < nx*order+1; i++)
    {
        tx[i] = (double)i/order;
    }
    double* ty = new double[ny*order+1];
    for (int i = 0; i < ny*order+1; i++)
    {
        ty[i] = (double)i/order;
    }
    const int ix = (int)flag/(order+1);
    const int iy = flag%(order+1);
    for (int j = 0; j < order+1; j++)
    {
        if (j != ix)
        {
            for (int k = 0; k < gs*gs*nx*ny; k++)
            {
                poly_term[k] *= (x_term[k] - tx[j])/(tx[ix] - tx[j]);
            }
        }
        if (j != iy)
        {
            for (int k = 0; k < gs*gs*nx*ny; k++)
            {
                poly_term[k] *= (y_term[k] - ty[j])/(ty[iy] - ty[j]);
            }
        }
    }
    

    delete [] pol_x1;
    delete [] pol_y1;
    delete [] x_term;
    delete [] y_term;
    delete [] tx;
    delete [] ty;

}