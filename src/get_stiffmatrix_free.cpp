#include "head.h"
//#include <iostream>
void get_stiffmatrix_free(int nnz, int* c1, int* c2, double* M, int* c1_free, int* c2_free, 
    double* M_free, int* diri_set, int& nnz_free)
{
    int ind = 0;
    for(int i = 0; i < nnz; i++)
    {
        if((diri_set[c1[i]-1] == -1) || (diri_set[c2[i]-1] == -1)) // 是Dirichlet点
        {
            continue;
        }
        else 
        {
            c1_free[ind] = c1[i] - diri_set[c1[i]-1];
            c2_free[ind] = c2[i] - diri_set[c2[i]-1];
            M_free[ind] = M[i];
            ind++;
        }
    }
    nnz_free = ind;
}