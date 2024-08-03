// get_b.cpp
#include <mkl.h>
#include "head.h"
#include <mkl_spblas.h>
using namespace std;
void get_b(const double x0, const double y0, const double xl, const double yl, 
    const int nx, const int ny, const int order, const double* inpol, 
    const double* coeff, const double* f, double* b)
{
    const int gs = order + 1;

    double hx = xl - x0;
    double hy = yl - y0;
    double* pol_x = new double[gs];
    double* pol_y = new double[gs];
    for (int i = 0; i < gs; i++) {
        pol_x[i] = inpol[i];
        pol_y[i] = inpol[i];
    }
        // 利用高斯插值点计算网格内的所需高斯积分插值点
    add_constant_to_vector(pol_x, gs, 1.0);
    add_constant_to_vector(pol_y, gs, 1.0);
    cblas_dscal(gs, 1.0/2*hx/nx, (double *)pol_x, 1);
    cblas_dscal(gs, 1.0/2*hy/ny, (double *)pol_y, 1);
    
    double* W_vector = new double [gs*gs]; 
    for (int i = 0; i < gs*gs; i++) {
        W_vector[i] = 0.0;
    }
    // 执行列向量 x 乘以行向量 y 的操作 A = alpha * x * y^T
    cblas_dger(CblasRowMajor, gs, gs, 1, coeff, 1, coeff, 1, 
        (double *)W_vector, gs);

    int* s = new int[gs*gs*nx*ny];
    for(int i = 0; i < nx*ny; i++) {
        for(int j = 0; j < gs*gs; j++) {
            s[i*gs*gs+j] = (i/ny)*order*(ny*order+1) + (i%ny)*order + 1 + (j/gs)*(ny*order+1) + j%gs;
        }
    }
    // 打印s
    // print_vector("s", s, gs*gs*nx*ny);
    double* t_x = new double[gs];
    double* t_y = new double[gs];
    for(int i = 0; i < gs; i++) {
        t_x[i] = hx/nx/order*i;
        t_y[i] = hy/ny/order*i;
    }
    double* eta_x = new double[(gs)*(gs)];
    double* eta_y = new double[(gs)*(gs)];
    for (int i = 0; i < gs; i++) {
        for (int j = 0; j < gs; j++) {
            eta_x[i*gs+j] = 1;
            eta_y[i*gs+j] = 1;
        }
    }// 初始化 eta_x 和 eta_y
    for(int i = 0; i < gs; i++) {
        for(int j = 0; j < gs; j++) {
            if(j != i) {
                for(int k = 0; k < gs; k++) {
                    eta_x[k*gs+i] *= (pol_x[k] - t_x[j])/(t_x[i] - t_x[j]);
                    eta_y[k*gs+i] *= (pol_y[k] - t_y[j])/(t_y[i] - t_y[j]);
                }
            }
        }
    }
    double eta[gs*gs][gs*gs];
    
    for(int i = 0; i < gs; i++) {
        for(int j = 0; j < gs; j++) {
            double* vecy = new double[gs];
            extract_column_as_vector((double*)eta_y, gs, gs, j, vecy);
            double* vecx = new double[gs];
            extract_column_as_vector((double*)eta_x, gs, gs, i, vecx);
            double* tmp1_0 = new double[gs*gs];
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, gs, gs, 1, 
                1.0, (double *)vecy, gs, (double *)vecx, gs, 0.0, tmp1_0, gs);
            for(int k = 0; k < gs*gs; k++) {
                eta[k][i*gs+j] = tmp1_0[k];
            }
            delete [] vecy;
            delete [] vecx;
            delete [] tmp1_0;
            
        }
    }
    double tmp3[gs*gs][gs*gs];
    vector_matrix_row_multiply((double*)W_vector, (double*)eta,
        (double*)tmp3, gs*gs, gs*gs);
    
    cblas_dscal(gs*gs*gs*gs, 1.0/4/nx/ny*hx*hy, (double*)tmp3, 1);
    // 定义B
    double* B = new double [gs*gs*nx*ny];
    for(int i = 0; i < nx*ny; i++)
    {
        const double* f_local = f + gs*gs*i;
        double tmpB[gs*gs][gs*gs];
        vector_matrix_row_multiply((double*)f_local, (double*)tmp3,
        (double*)tmpB, gs*gs, gs*gs);
        double* sumB = new double[gs*gs];
        matrix_column_sum((double*)tmpB, gs*gs, gs*gs, sumB);
        for(int j = 0; j < gs*gs; j++)
        {
            B[i*gs*gs+j] = sumB[j];
        }
        delete [] sumB;

    }
    for(int i = 0; i < (nx*order+1)*(ny*order+1); i++)
    {
        b[i] = 0;
    }
    for(int i = 0; i < gs*gs*nx*ny; i++)
    {
        b[s[i]-1] += B[i];
    }
    delete [] s;
    delete [] pol_x;
    delete [] pol_y;
    delete [] W_vector;
    delete [] t_x;
    delete [] t_y;
    delete [] B;
    delete [] eta_x;
    delete [] eta_y;
    
    
    
}