#include "head.h"
#include <iostream>
#include <cmath>
using namespace std;
void get_poly_bound(const int order, const int nx, const int ny, double *pbd, const int flag)
{
    if (flag < 0 || flag >= (order + 1) * (order + 1))    // 从0开始计数
    {
        cout << "flag is out of range" << endl;
    }
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
    double *tx = new double[order + 1];
    for (int i = 0; i < order + 1; i++)
    {
        tx[i] = (double)i / order;
    }
    double *ty = new double[order + 1];
    for (int i = 0; i < order + 1; i++)
    {
        ty[i] = (double)i / order;
    }
    const int ix = floor(flag / (order + 1));
    const int iy = flag % (order + 1);
    for (int i = 0; i < sizebd; i++)
    {
        pbd[i] = 1.0;
    }
    for (int j = 0; j < order + 1; j++)
    {
        if (j != ix)
        {
            for (int k = 0; k < sizebd; k++)
            {
                pbd[k] *= (bdx[k] - tx[j]) / (tx[ix] - tx[j]);
            }
        }
        if (j != iy)
        {
            for (int k = 0; k < sizebd; k++)
            {
                pbd[k] *= (bdy[k] - ty[j]) / (ty[iy] - ty[j]);
            }
        }
    }
    delete[] tx;
    delete[] ty;
    delete[] bdx;
    delete[] bdy;
}