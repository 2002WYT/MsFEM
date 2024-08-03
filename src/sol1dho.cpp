#include "head.h"
#include <iostream>
#include <mkl.h>
using namespace std;
void sol1dho(const double x0, const double xl, const double* coeff, const double* inpol,
    const int n, const int order, const double* bound_val, const double* k, const double* f, 
    double* u)
{
    const double h = xl - x0;
    const int gs = order + 1;
    double* pol = new double[gs];
    for (int i = 0; i < gs; i++) {
        pol[i] = (inpol[i] + 1.0)/2.0*h/n;
    }
    double* t = new double[gs];
    for (int i = 0; i < gs; i++) {
        t[i] = i * h / n / order;
    }
    double* eta = new double[gs*gs];
    double* d_eta = new double[gs*gs];
    for (int i = 0; i < gs*gs; i++) {
        eta[i] = 1.0;
        d_eta[i] = 0.0;
    }
    for (int i = 0; i < gs; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            if (j!=i)
            {
                for (int k = 0; k < gs; k++)
                {
                    eta[i*gs+k] *= (pol[k] - t[j]) / (t[i] - t[j]); // eta列存储
                }
                double* tmp = (double*)mkl_malloc(gs*sizeof(double), 64);
                for (int k = 0; k < gs; k++)
                {
                    tmp[k] = 1.0/(t[i] - t[j]);}
                for (int p = 0; p < gs; p++)
                {
                if (p!=i && p!=j)
                    {
                        for (int k = 0; k < gs; k++)
                        {
                            tmp[k] *= (pol[k] - t[p]) / (t[i] - t[p]);
                        }
                    }
                }
                for (int k = 0; k < gs; k++)
                {
                    d_eta[i*gs+k] += tmp[k];
                }
                mkl_free(tmp);
            }
        }
    }
    int* s = new int[n*gs];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            s[i*gs+j] = i * order + j + 1;
        }
    }
    double* J = new double[gs*(gs*gs)];
    for (int i = 0; i < gs; i++)
    {
        for (int j = 0; j < gs*gs; j++)
        {
            J[j*gs + i] = d_eta[j/gs*gs + i] * d_eta[j%gs*gs+i]; // 列存储
        }
    }
    double* tmp2c = new double[gs];
    for (int i = 0; i < gs; i++) tmp2c[i] = coeff[i];
    cblas_dscal(gs, 1.0/2*h/n, tmp2c, 1);
    double* tmpe = new double[gs*gs];
    for (int i = 0; i < gs*gs; i++)
    {
        tmpe[i] = eta[i] * tmp2c[i%gs];
    }
    double* tmpJ = new double[gs*(gs*gs)];
    for (int i = 0; i < gs*(gs*gs); i++)
    {
        tmpJ[i] = J[i] * tmp2c[i%gs];
    }
    double* M = new double[gs*gs*n];
    double* B = new double[gs*n];
    for (int i = 0; i < n; i++)
    {
        double* tmpJ2 = new double[gs*(gs*gs)];
        double* tmpe2 = new double[gs*gs];
        double* k_local = new double[gs];
        double* f_local = new double[gs];
        for (int j = 0; j < gs; j++)
        {
            k_local[j] = k[i*gs+j];
            f_local[j] = f[i*gs+j];
        }
        for (int j = 0; j < gs*(gs*gs); j++)
        {
            tmpJ2[j] = tmpJ[j] * k_local[j%gs];
        }
        for (int j = 0; j < gs*gs; j++)
        {
            tmpe2[j] = tmpe[j] * f_local[j%gs];
        }
        double* sumM = new double[gs*gs];
        double* sumB = new double[gs];
        matrixcol_column_sum(tmpJ2, gs, gs*gs, sumM);
        matrixcol_column_sum(tmpe2, gs, gs, sumB);
        for (int j = 0; j < gs*gs; j++)
        {
            M[i*gs*gs+j] = sumM[j];
        }
        for (int j = 0; j < gs; j++)
        {
            B[i*gs+j] = sumB[j];
        }
        
        delete [] k_local;
        delete [] f_local;
        delete [] sumM;
        delete [] sumB;
        delete [] tmpJ2;
        delete [] tmpe2;
    }
    int* c1 = new int [n*gs*gs];
    int* c2 = new int [n*gs*gs];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < gs*gs; j++)
        {
            c1[i*gs*gs+j] = s[i*gs+j/gs];
            c2[i*gs*gs+j] = s[i*gs+j%gs];
        }
    }
    double* b = new double [n*order+1];
    for(int i = 0; i < n*order+1; i++)
    {
        b[i] = 0;
    }
    for(int i = 0; i < gs*n; i++)
    {
        b[s[i]-1] += B[i];
    }
    b[0] = bound_val[0];
    b[n*order] = bound_val[1];
    delete [] s;
    // 调整矩阵A
    for (int i = 0; i < n*gs*gs; i++)
    {
        if (c1[i] == 1 || c1[i] == n*order+1)
        {
            if (c2[i] == c1[i]){
                M[i] = 1.0;
            }
            else M[i] = 0.0;
        }
    }
    const int nnz = n*gs*gs;
    const int sizeA = n*order+1;

    sparse_matrix_t A;
    mkl_sparse_d_create_coo(&A, SPARSE_INDEX_BASE_ONE,
        sizeA, sizeA, nnz, c1, c2, M);
    sparse_matrix_t csr_A;
    mkl_sparse_convert_csr(A, SPARSE_OPERATION_NON_TRANSPOSE, &csr_A);
    sparse_index_base_t indexing;
    int rows_csr, cols_csr;
    int* row_ptr;
    int* col_ind;
    int* row_end;
    double* csr_val;
    mkl_sparse_d_export_csr(csr_A, &indexing, &rows_csr, &cols_csr, 
        &row_ptr, &row_end, &col_ind, &csr_val);
    
    delete [] pol;
    delete [] t;
    delete [] eta;
    delete [] d_eta;
    delete [] c1;
    delete [] c2;
    delete [] J;
    delete [] tmp2c;
    delete [] tmpe;
    delete [] tmpJ;
    delete [] M;
    delete [] B;
    
    // 定义PARDISO求解器的参数
    void *pt[64] = {0};
    // 控制参数数组
    MKL_INT iparm[64] = {0};

    // PARDISO求解器的其他参数
    MKL_INT maxfct = 1;  // 求解不同矩阵时设置为1
    MKL_INT mnum = 1;    // 求解单个矩阵时设置为1
    MKL_INT mtype = 1;  // 实数对称矩阵
    MKL_INT phase = 13;  // 求解阶段，初始化为13
    MKL_INT nrhs = 1;    // 右手边向量数量
    MKL_INT error = 0;   // 错误标志
    MKL_INT msglvl = 0;  // 输出级别

    // 初始化iparm数组，设置PARDISO的默认参数
    iparm[0] = 1;  // 使用默认值
    iparm[1] = 2;  // 并行处理
    iparm[5] = 0;  // 
    iparm[7] = 2;  // 迭代步骤
    iparm[9] = 8; // 对称矩阵
    iparm[10] = 1; // 稳定的默认设置
    // 调用PARDISO求解器!!!
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &sizeA, 
        csr_val, row_ptr, col_ind, NULL, &nrhs, iparm, 
        &msglvl, b, u, &error);

    // 检查是否有错误
    if (error != 0) {
        std::cerr << "Error during solution: " << error << std::endl;
    }
    
    // 打印解向量
    // print_vector("u", u, sizeA);
    phase = -1;  // 释放内存阶段
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &sizeA, nullptr, 
        row_ptr, col_ind, nullptr, &nrhs, iparm, &msglvl, 
        nullptr, nullptr, &error);

    
    delete [] b;

    mkl_sparse_destroy(csr_A);
    mkl_sparse_destroy(A);
}