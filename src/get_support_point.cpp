#include "head.h"
#include <mkl.h>
void get_support_point(const int nx, const int ny, const int order, int* suppind)
{
    const int gs = order+1;
    for (int i = 0; i < gs; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            suppind[i*gs+j] = i*(ny*order+1)*nx+j*ny;
        }
    }
}