// 一次性求gs*gs个振荡边界条件，osbd大小为gs*gs*(2*ny*order+2*nx*order)
#include "head.h"
// #include <iostream>
#include <mkl.h>
#include <cmath>
using namespace std;
void get_oscillatory_bound_2(const double x0, const double xl, const double y0, const double yl,
                           const double *k_column1, const double *k_column2, const double *k_row1, const double *k_row2,
                           const int nx, const int ny, const int order, const double *coeff, const double *inpol,
                           double *osbd)
{
    const int gs = order + 1;

    const double Lx    = xl - x0;
    const double Ly    = yl - y0;
    double      *pol_x = new double[gs];
    double      *pol_y = new double[gs];
    for (int i = 0; i < gs; i++)
    {
        pol_x[i] = (inpol[i] + 1.0) / 2 * Lx / nx;
        pol_y[i] = (inpol[i] + 1.0) / 2 * Ly / ny;
    }
    // 利用高斯插值点计算网格内的所需高斯积分插值点
    double *xpnt = new double[gs * nx];
    double *ypnt = new double[gs * ny];
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            xpnt[i * gs + j] = i * Lx / nx + pol_x[j];
        }
    }
    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            ypnt[i * gs + j] = i * Ly / ny + pol_y[j];
        }
    }
    double *tx = new double[gs];
    for (int i = 0; i < gs; i++)
    {
        tx[i] = (double)i * Lx / order;
    }
    double *ty = new double[gs];
    for (int i = 0; i < gs; i++)
    {
        ty[i] = (double)i * Ly / order;
    }
    double *etax = new double[(gs * nx) * gs];
    double *etay = new double[(gs * ny) * gs];
    // 一次循环初始化为全1
    for (int i = 0; i < gs * nx * gs; i++)
    {
        etax[i] = 1.0;
    }
    for (int i = 0; i < gs * ny * gs; i++)
    {
        etay[i] = 1.0;
    }

    for (int i = 0; i < gs; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            if (j != i)
            {
                for (int k = 0; k < gs * nx; k++)
                {
                    etax[i * gs * nx + k] *= (xpnt[k] - tx[j]) / (tx[i] - tx[j]);
                }
                for (int k = 0; k < gs * ny; k++)
                {
                    etay[i * gs * ny + k] *= (ypnt[k] - ty[j]) / (ty[i] - ty[j]);
                }
            }
        }
    }
    double *gamma1 = new double[(ny * order + 1) * gs];
    // 单元左边界的一维振荡边界条件，也可看作列优先存储的矩阵。
    double *gamma2 = new double[(nx * order + 1) * gs];    // 单元下边界的一维振荡边界条件
    double *gamma3 = new double[(nx * order + 1) * gs];    // 单元上边界的一维振荡边界条件
    double *gamma4 = new double[(ny * order + 1) * gs];    // 单元右边界的一维振荡边界条件
    for (int i = 0; i < (ny * order + 1) * gs; i++)
    {
        gamma1[i] = 0.0;
        gamma4[i] = 0.0;
    }
    for (int i = 0; i < (nx * order + 1) * gs; i++)
    {
        gamma2[i] = 0.0;
        gamma3[i] = 0.0;
    }
    double  db1[2]  = { 1.0, 0.0 };
    double  db2[2]  = { 0.0, 1.0 };
    double  db3[2]  = { 0.0, 0.0 };
    double *zerocol = new double[gs * ny];
    double *zerorow = new double[gs * nx];
    for (int i = 0; i < gs * ny; i++)
        zerocol[i] = 0.0;
    for (int i = 0; i < gs * nx; i++)
        zerorow[i] = 0.0;

    for (int i = 0; i < gs; i++)
    {
        int indexy = i * (ny * order + 1);
        int indexx = i * (nx * order + 1);
        if (i == 0)
        {
            sol1dho(y0, yl, coeff, inpol, ny, order, db1, k_column1,
                    zerocol, gamma1 + indexy);
            sol1dho(x0, xl, coeff, inpol, nx, order, db1, k_row1,
                    zerorow, gamma2 + indexx);
            sol1dho(x0, xl, coeff, inpol, nx, order, db1, k_row2,
                    zerorow, gamma3 + indexx);
            sol1dho(y0, yl, coeff, inpol, ny, order, db1, k_column2,
                    zerocol, gamma4 + indexy);
        }
        else if (i == gs - 1)
        {
            sol1dho(y0, yl, coeff, inpol, ny, order, db2, k_column1,
                    zerocol, gamma1 + indexy);
            sol1dho(x0, xl, coeff, inpol, nx, order, db2, k_row1,
                    zerorow, gamma2 + indexx);
            sol1dho(x0, xl, coeff, inpol, nx, order, db2, k_row2,
                    zerorow, gamma3 + indexx);
            sol1dho(y0, yl, coeff, inpol, ny, order, db2, k_column2,
                    zerocol, gamma4 + indexy);
        }
        else
        {
            sol1dho(y0, yl, coeff, inpol, ny, order, db3, k_column1,
                    etay + i * gs * ny, gamma1 + indexy);
            sol1dho(x0, xl, coeff, inpol, nx, order, db3, k_row1,
                    etax + i * gs * nx, gamma2 + indexx);
            sol1dho(x0, xl, coeff, inpol, nx, order, db3, k_row2,
                    etax + i * gs * nx, gamma3 + indexx);
            sol1dho(y0, yl, coeff, inpol, ny, order, db3, k_column2,
                    etay + i * gs * ny, gamma4 + indexy);
        }
    }
    double *gsp1 = new double[gs * gs];
    double *gsp2 = new double[gs * gs];
    double *gsp3 = new double[gs * gs];
    double *gsp4 = new double[gs * gs];    // 都是各一维基函数在支撑点的值列优先存储。
    for (int i = 0; i < gs; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            gsp1[i * gs + j] = gamma1[i * (ny * order + 1) + j * ny];
            gsp2[i * gs + j] = gamma2[i * (nx * order + 1) + j * nx];
            gsp3[i * gs + j] = gamma3[i * (nx * order + 1) + j * nx];
            gsp4[i * gs + j] = gamma4[i * (ny * order + 1) + j * ny];
        }
    }
    int *ipiv1 = (int *)mkl_malloc(gs * sizeof(int), 64);
    // 调用dgetrf进行LU分解
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, gs, gs, gsp1, gs, ipiv1);
    // 调用dgetri计算逆矩阵
    LAPACKE_dgetri(LAPACK_COL_MAJOR, gs, gsp1, gs, ipiv1);

    int *ipiv2 = (int *)mkl_malloc(gs * sizeof(int), 64);
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, gs, gs, gsp2, gs, ipiv2);
    LAPACKE_dgetri(LAPACK_COL_MAJOR, gs, gsp2, gs, ipiv2);
    int *ipiv3 = (int *)mkl_malloc(gs * sizeof(int), 64);
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, gs, gs, gsp3, gs, ipiv3);
    LAPACKE_dgetri(LAPACK_COL_MAJOR, gs, gsp3, gs, ipiv3);
    int *ipiv4 = (int *)mkl_malloc(gs * sizeof(int), 64);
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, gs, gs, gsp4, gs, ipiv4);
    LAPACKE_dgetri(LAPACK_COL_MAJOR, gs, gsp4, gs, ipiv4);
    // 这里求逆比MATLAB精度高很多
    mkl_free(ipiv1);
    mkl_free(ipiv2);
    mkl_free(ipiv3);
    mkl_free(ipiv4);
    // print_matrix_col("gsp1", gsp1, gs, gs);
    double *e1 = new double[(ny * order + 1) * gs];
    double *e2 = new double[(nx * order + 1) * gs];
    double *e3 = new double[(nx * order + 1) * gs];
    double *e4 = new double[(ny * order + 1) * gs];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ny * order + 1,
                gs, gs, 1.0, gamma1,
                ny * order + 1, gsp1, gs, 0.0, e1, ny * order + 1);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nx * order + 1,
                gs, gs, 1.0, gamma2,
                nx * order + 1, gsp2, gs, 0.0, e2, nx * order + 1);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nx * order + 1,
                gs, gs, 1.0, gamma3,
                nx * order + 1, gsp3, gs, 0.0, e3, nx * order + 1);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ny * order + 1,
                gs, gs, 1.0, gamma4,
                ny * order + 1, gsp4, gs, 0.0, e4, ny * order + 1);

    delete[] pol_x;
    delete[] pol_y;
    delete[] xpnt;
    delete[] ypnt;
    delete[] tx;
    delete[] ty;
    delete[] etax;
    delete[] etay;
    delete[] gamma1;
    delete[] gamma2;
    delete[] gamma3;
    delete[] gamma4;
    delete[] gsp1;
    delete[] gsp2;
    delete[] gsp3;
    delete[] gsp4;
    delete[] zerocol;
    delete[] zerorow;
    int sizebd = 2 * ny * order + 2 * nx * order;

    for (int flag = 0; flag < gs * gs; flag++)
    {
        int flindex = flag * sizebd;
        if (flag == 0)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = e1[i];
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                int index = i - (ny * order + 1) + 2;
                if (index % 2 == 0)
                    osbd[i + flindex] = e2[index / 2];
                else
                    osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                osbd[i + flindex] = 0.0;
            }
        }
        else if (flag == gs - 1)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = e1[i + order * (ny * order + 1)];
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                int index = i - (ny * order + 1) + 1;
                if (index % 2 == 0)
                    osbd[i + flindex] = e3[index / 2];
                else
                    osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                osbd[i + flindex] = 0.0;
            }
        }
        else if (flag == gs * gs - gs)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                int index = i - (ny * order + 1) + 2;
                if (index % 2 == 0)
                    osbd[i + flindex] = e2[order * (nx * order + 1) + index / 2];
                else
                    osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                int index = i - (ny * order + 2 * nx * order - 1);
                osbd[i + flindex]   = e4[index];
            }
        }
        else if (flag == gs * gs - 1)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                int index = i - (ny * order + 1) + 1;
                if (index % 2 == 0)
                    osbd[i + flindex] = e3[order * (nx * order + 1) + index / 2];
                else
                    osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                int index = i - (ny * order + 2 * nx * order - 1);
                osbd[i + flindex]   = e4[order * (ny * order + 1) + index];
            }
        }
        else if (flag > 0 && flag < gs - 1)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = e1[flag * (ny * order + 1) + i];
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                osbd[i + flindex] = 0.0;
            }
        }
        else if (flag > gs - 1 && flag < gs * gs - gs - 1 && flag % gs == 0)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                int index = i - (ny * order + 1) + 2;
                int bind  = (int)flag / gs;
                if (index % 2 == 0)
                    osbd[i + flindex] = e2[bind * (nx * order + 1) + index / 2];
                else
                    osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                osbd[i + flindex] = 0.0;
            }
        }
        else if (flag > gs && flag < gs * gs - 1 && flag % gs == order)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                int index = i - (ny * order + 1) + 1;
                int bind  = (int)flag / gs;
                if (index % 2 == 0)
                    osbd[i + flindex] = e3[bind * (nx * order + 1) + index / 2];
                else
                    osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                osbd[i + flindex] = 0.0;
            }
        }
        else if (flag > gs * gs - gs && flag < gs * gs - 1)
        {
            for (int i = 0; i < ny * order + 1; i++)
            {
                osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 1; i < ny * order + 2 * nx * order - 1; i++)
            {
                osbd[i + flindex] = 0.0;
            }
            for (int i = ny * order + 2 * nx * order - 1; i < sizebd; i++)
            {
                int index = i - (ny * order + 2 * nx * order - 1);
                osbd[i + flindex]   = e4[(flag - gs * order) * (ny * order + 1) + index];
            }
        }
        else // 内部
        {
            for (int i = 0; i < sizebd; i++)
            {
                osbd[i + flindex] = 0.0;
            }
        }
    }

    delete[] e1;
    delete[] e2;
    delete[] e3;
    delete[] e4;
}
