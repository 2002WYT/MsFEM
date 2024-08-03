
#include <mkl.h>
#include "head.h"
#include <mkl_spblas.h>
#include <iostream>
using namespace std;
void get_stiffmatrix(const double x0, const double y0, const double xl, const double yl, 
    const int nx, const int ny, const int order, const double* inpol, 
    const double* coeff, const double* k, int* c1, int* c2, double* M)
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
        W_vector[i] = 0;
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
    double* dex = new double[(gs)*(gs)];
    double* dey = new double[(gs)*(gs)];
    for (int i = 0; i < gs; i++) {
        for (int j = 0; j < gs; j++) {
            dex[i*gs+j] = 0;
            dey[i*gs+j] = 0;
        }
    }// 初始化 dex 和 dey
    for(int i = 0; i < gs; i++) {
        for(int j = 0; j < gs; j++) {
            if(j != i) {
                for(int k = 0; k < gs; k++) {
                    eta_x[k*gs+i] *= (pol_x[k] - t_x[j])/(t_x[i] - t_x[j]);
                    eta_y[k*gs+i] *= (pol_y[k] - t_y[j])/(t_y[i] - t_y[j]);
                }
                double* tmp_x = new double[gs];
                double* tmp_y = new double[gs];
                for(int k = 0; k < gs; k++) {
                    tmp_x[k] = 1.0 / (t_x[i] - t_x[j]);
                    tmp_y[k] = 1.0 / (t_y[i] - t_y[j]);
                }
                for(int p = 0; p < gs; p++) {
                    if(p != i && p != j) {
                        for(int k = 0; k < gs; k++) {
                            tmp_x[k] *= (pol_x[k] - t_x[p])/(t_x[i] - t_x[p]);
                            tmp_y[k] *= (pol_y[k] - t_y[p])/(t_y[i] - t_y[p]);
                        }
                    }
                }

                for(int k = 0; k < gs; k++) {
                    dex[k*gs+i] += tmp_x[k];
                    dey[k*gs+i] += tmp_y[k];
                }
                delete [] tmp_x;
                delete [] tmp_y;
            }
        }
    }
    double d_eta_x[gs*gs][gs*gs], d_eta_y[gs*gs][gs*gs];
    for(int i = 0; i < gs; i++) {
        for(int j = 0; j < gs; j++) {
            double* vecy = new double[gs];
            extract_column_as_vector((double*)eta_y, gs, gs, j, vecy);
            double* vecx = new double[gs];
            extract_column_as_vector((double*)eta_x, gs, gs, i, vecx);
           
            double* vecd1 = new double[gs];
            extract_column_as_vector((double*)dex, gs, gs, i, vecd1);// 提取dex的第i列
            double* vecd2 = new double[gs];
            extract_column_as_vector((double*)dey, gs, gs, j, vecd2);// 提取dex的第j列
            double* tmp1_1 = new double[gs*gs];
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, gs, gs, 1, 
                1.0, (double *)vecy, gs, (double *)vecd1, gs, 0.0, tmp1_1, gs);
            for(int k = 0; k < gs*gs; k++) {
                d_eta_x[k][i*gs+j] = tmp1_1[k];
            }
            double* tmp1_2 = new double[gs*gs];
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, gs, gs, 1, 
                1.0, (double *)vecd2, gs, (double *)vecx, gs, 0.0, tmp1_2, gs);
            for(int k = 0; k < gs*gs; k++) {
                d_eta_y[k][i*gs+j] = tmp1_2[k];
            }
            delete [] tmp1_1;
            delete [] tmp1_2;
            delete [] vecy;
            delete [] vecx;
            delete [] vecd1;
            delete [] vecd2;
        }
    }

    double J[gs*gs][gs*gs*gs*gs];
    for(int i = 0; i < gs*gs*gs*gs; i++) {
        for(int j = 0; j < gs*gs; j++) {
            J[j][i] = d_eta_x[j][i/(gs*gs)]*d_eta_x[j][i%(gs*gs)] + d_eta_y[j][i/(gs*gs)]*d_eta_y[j][i%(gs*gs)];
        }
    }
    double tmp2[gs*gs][gs*gs*gs*gs];
    vector_matrix_row_multiply((double*)W_vector, (double*)J,
        (double*)tmp2, gs*gs, gs*gs*gs*gs);
    cblas_dscal(gs*gs*gs*gs*gs*gs, 1.0/4/nx/ny*hx*hy, (double*)tmp2, 1);
    
    // 定义M和B
    for(int i = 0; i < nx*ny; i++)
    {
        const double* k_local = k + gs*gs*i;
        
        double* tmpM = new double [(gs*gs)*(gs*gs*gs*gs)];
        
        vector_matrix_row_multiply((double*)k_local, (double*)tmp2,
        (double*)tmpM, gs*gs, gs*gs*gs*gs); 
        double* sumM = new double[gs*gs*gs*gs];
        for (int j = 0; j < gs*gs*gs*gs; j++) {
            sumM[j] = 0;
        }
        
        matrix_column_sum(tmpM, gs*gs, gs*gs*gs*gs, sumM);
        for(int j = 0; j < gs*gs*gs*gs; j++)
        {
            M[i*gs*gs*gs*gs+j] = sumM[j];
        }
        
        delete [] tmpM;
        delete [] sumM;
    }
    // 将s每行复制gs份拼凑到一起（张量积）
    for(int i = 0; i < nx*ny; i++)
    {
        for(int j = 0; j < gs*gs*gs*gs; j++)
        {
            c1[i*gs*gs*gs*gs+j] = s[i*gs*gs+j/(gs*gs)];
            c2[i*gs*gs*gs*gs+j] = s[i*gs*gs+j%(gs*gs)];
        }
    }
    delete [] s;
    delete [] pol_x;
    delete [] pol_y;
    delete [] W_vector;
    delete [] t_x;
    delete [] t_y;
    delete [] eta_x;
    delete [] eta_y;
    delete [] dex;
    delete [] dey;

    
    
}