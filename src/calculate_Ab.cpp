// subfunction
#include <mkl.h>
#include "head.h"
#include <mkl_spblas.h>
// #include <chrono>
// #include <iostream>
using namespace std;
void calculate_Ab(const int nx, const int ny, const int order, const double* bound_val, int nnz, 
    int* c1, int* c2, double* M, double* b, double* u)
{
    const int sizeA = (nx*order+1)*(ny*order+1);
    sparse_matrix_t A;
    mkl_sparse_d_create_coo(&A, SPARSE_INDEX_BASE_ONE,
        sizeA, sizeA, nnz, c1, c2, M);
    sparse_matrix_t csr_A;
    mkl_sparse_convert_csr(A, SPARSE_OPERATION_NON_TRANSPOSE, &csr_A);
    const int sizefree = sizeA - 2*ny*order - 2*nx*order;
    sparse_index_base_t indexing;
    int rows_csr, cols_csr;
    int* row_ptr; 
    int* col_ind;
    int* row_end;
    double* csr_val;
    mkl_sparse_d_export_csr(csr_A, &indexing, &rows_csr, &cols_csr, 
        &row_ptr, &row_end, &col_ind, &csr_val);
    
    int* diri_set = new int [sizeA]; // 存储Dirichlet点的位置
    double* diri_val = new double [sizeA]; // 存储Dirichlet点的值
    int cnt = 0;
    for(int i = 0; i < sizeA; i++)
    {
        if(i < ny*order+1 || i % (ny*order+1) == 0 || i % (ny*order+1) == ny*order ||
            i > (nx*order+1)*(ny*order+1)-ny*order-1)
        {
            diri_set[i] = -1; // 是Dirichlet点
            diri_val[i] = bound_val[cnt];
            cnt++;
        }
        else
        {
            diri_set[i] = cnt;
            diri_val[i] = 0; // 不是Dirichlet点
        }
    }
    int* c1_free = new int [nnz];
    int* c2_free = new int [nnz];
    double* M_free = new double [nnz];
    int nnz_free = 0;
    get_stiffmatrix_free(nnz, c1, c2, M, c1_free, c2_free, M_free, diri_set, nnz_free);
    sparse_matrix_t A_free;
    mkl_sparse_d_create_coo(&A_free, SPARSE_INDEX_BASE_ONE,
        sizeA, sizeA, nnz_free, c1_free, c2_free, M_free);
    sparse_matrix_t csr_A_free;
    mkl_sparse_convert_csr(A_free, SPARSE_OPERATION_NON_TRANSPOSE, &csr_A_free);
    
    sparse_index_base_t indexing_free;
    int rows_csr_free, cols_csr_free;
    int *row_ptr_free;
    int* col_ind_free;
    int* row_end_free;
    double *csr_val_free;
    mkl_sparse_d_export_csr(csr_A_free, &indexing_free, &rows_csr_free, 
        &cols_csr_free, &row_ptr_free, &row_end_free, 
        &col_ind_free, &csr_val_free);

    for(int i = 0; i < sizeA; i++)
    {
        if(diri_set[i] != -1)
        {
            u[i] = 0; // 非Dirichlet边界值
        }
        else
        {
            u[i] = diri_val[i]; // Dirichlet边界值
        }
    }
    // 执行b = b - A*u操作
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, csr_A, 
        descr, u, 1.0, b);
    double* b_free = new double [sizefree];
    double* u_free = new double [sizefree]; // 存储自由点的解
    get_rightterm_free(sizeA, b, diri_set, b_free); // 得到自由点的右端项
    for(int i = 0; i < sizefree; i++)
    {
        u_free[i] = 0;
    }
    // 打印diri_set
    //print_vector("diri_set", diri_set, sizeA);
    // 在求解前把无用的内存释放掉
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
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &sizefree, 
        csr_val_free, row_ptr_free, col_ind_free, NULL, &nrhs, iparm, 
        &msglvl, b_free, u_free, &error);
    cnt = 0;
    for(int i = 0; i < sizeA; i++)
    {
        if(diri_set[i] == -1) // 是Dirichlet点
        {
            continue;
        }
        else
        {
            u[i] = u_free[cnt];
            cnt++;
        }
    }
    // 打印解向量
    // print_vector("u", u, sizeA);
    phase = -1;  // 释放内存阶段
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &sizefree, nullptr, 
        row_ptr_free, col_ind_free, nullptr, &nrhs, iparm, &msglvl, 
        nullptr, nullptr, &error);
    // 释放内存
    mkl_sparse_destroy(A_free);
    mkl_sparse_destroy(csr_A_free);
    mkl_sparse_destroy(A);
    mkl_sparse_destroy(csr_A);
    delete [] diri_val;
    delete [] diri_set;
    delete [] u_free;
    //delete [] c1;
    //delete [] c2;
    //delete [] M;
    delete [] c1_free;
    delete [] c2_free;
    delete [] M_free;
    delete [] b_free;

    
    

}