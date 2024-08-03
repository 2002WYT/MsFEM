#include "head.h"
void get_s_fine(const int Nx, const int Ny, const int nx, const int ny, const int order, int* s_fine)
{
    for(int i = 0; i < Nx*Ny; i++)
    {
        for(int j = 0; j < (nx*order+1)*(ny*order+1); j++)
        {
            s_fine[i*(nx*order+1)*(ny*order+1)+j] = (i/Ny)*nx*order*(Ny*ny*order+1) + 
            (i%Ny)*ny*order + 1 + (j/(ny*order+1))*(Ny*ny*order+1) + j%(ny*order+1);
        }
    }
}
