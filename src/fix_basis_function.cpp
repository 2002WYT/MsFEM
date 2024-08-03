#include "head.h"
#include <mkl.h>
void fix_basis_function(const int nx, const int ny, const int order, const double* R_local, 
    double* R_fix)
{
    const int gs = order+1;
    double* R_supp = (double*)mkl_malloc((gs*gs)*(gs*gs)*sizeof(double), 64);
    int* suppind = (int *)mkl_malloc(gs*gs*sizeof(int), 64);
    get_support_point(nx, ny, order, suppind);
    for (int i = 0; i < gs*gs; i++)
    {
        for (int j = 0; j < gs*gs; j++)
        {
            R_supp[i*gs*gs+j] = R_local[i*(ny*order+1)*(nx*order+1)+suppind[j]];
        }
    }
    int* ipiv1 = (int *)mkl_malloc(gs*gs*sizeof(int), 64);
    // 调用dgetrf进行LU分解
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, gs*gs, gs*gs, R_supp, gs*gs, 
        ipiv1);
    // 调用dgetri计算逆矩阵
    LAPACKE_dgetri(LAPACK_COL_MAJOR, gs*gs, R_supp, gs*gs, ipiv1);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
        (ny*order+1)*(nx*order+1), gs*gs, gs*gs, 1.0, R_local, 
        (ny*order+1)*(nx*order+1), R_supp, gs*gs, 0.0, R_fix, 
        (ny*order+1)*(nx*order+1));





    mkl_free(R_supp);
    mkl_free(suppind);
    mkl_free(ipiv1);

}