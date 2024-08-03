// 一次性求gs*gs个多项式边界条件，pbd大小为sizebd*gs*gs
#include "head.h"
// #include <iostream>
#include <cmath>
using namespace std;
void get_poly_bound_2(const int order, const int nx, const int ny, double *pbd)
{
    const int gs     = order + 1;
    const int sizebd = 2 * nx * order + 2 * ny * order;
    double   *bdx    = new double[sizebd];
    double   *bdy    = new double[sizebd];
    for (int i = 0; i < sizebd; i++)
    {
        if (i < ny * order + 1)
        {
            bdx[i] = 0.0;
            bdy[i] = (double)i / (ny * order);
        }
        else if (i >= ny * order + 1 && i < ny * order + 2 * nx * order - 1)
        {
            const int tmp = (int)(i - ny * order + 1) / 2;
            bdx[i]        = (double)tmp / (nx * order);
            bdy[i]        = (i - ny * order - 1) % 2;
        }
        else if (i >= ny * order + 2 * nx * order - 1 && i < sizebd)
        {

            bdx[i] = 1.0;
            bdy[i] = (double)(i - ny * order - 2 * nx * order + 1) / (ny * order);
        }
    }
    double *tx = new double[gs];
    for (int i = 0; i < gs; i++)
    {
        tx[i] = (double)i / order;
    }
    double *ty = new double[gs];
    for (int i = 0; i < gs; i++)
    {
        ty[i] = (double)i / order;
    }
    //
    for (int i = 0; i < sizebd * gs * gs; i++)
    {
        pbd[i] = 1.0;
    }
    for (int flag = 0; flag < gs * gs; flag++)
    {
        const int ix      = floor(flag / gs);
        const int iy      = flag % gs;
        int       flindex = flag * sizebd;

        for (int j = 0; j < gs; j++)
        {
            if (j != ix)
            {
                for (int k = 0; k < sizebd; k++)
                {
                    pbd[k + flindex] *= (bdx[k] - tx[j]) / (tx[ix] - tx[j]);
                }
            }
            if (j != iy)
            {
                for (int k = 0; k < sizebd; k++)
                {
                    pbd[k + flindex] *= (bdy[k] - ty[j]) / (ty[iy] - ty[j]);
                }
            }
        }
    }
    delete[] tx;
    delete[] ty;
    delete[] bdx;
    delete[] bdy;
}