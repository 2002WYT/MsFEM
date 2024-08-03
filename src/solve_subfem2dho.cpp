// subfunction
#include <mkl.h>
#include "head.h"
#include <mkl_spblas.h>
using namespace std;
void solve_subfem2dho(const double x0, const double y0, const double xl, const double yl, 
    const int nx, const int ny, const int order, const double* bound_val, const double* inpol, 
    const double* coeff, const double* k, const double* f, double* u)
{
    double hx = xl - x0;
    double hy = yl - y0;
    const int gs = order + 1;
    const int sizeA = (nx*order+1)*(ny*order+1);
    const int sizefree = sizeA - 2*ny*order - 2*nx*order;
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
    const int nnz = gs*gs*gs*gs*nx*ny;
    
    double* W_vector = new double [gs*gs]; 
    for (int i = 0; i < gs*gs; i++) {
        W_vector[i] = 0;
    } 
    // 执行列向量 x 乘以行向量 y 的操作 A = alpha * x * y^T
    cblas_dger(CblasRowMajor, gs, gs, 1, coeff, 1, coeff, 1, 
        (double *)W_vector, gs); // W_vector = coeff * coeff^T

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
    }
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
                    tmp_x[k] = 1 / (t_x[i] - t_x[j]);
                    tmp_y[k] = 1 / (t_y[i] - t_y[j]);
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
    double eta[gs*gs][gs*gs], d_eta_x[gs*gs][gs*gs], d_eta_y[gs*gs][gs*gs];
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
            delete [] tmp1_0;
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
    double tmp3[gs*gs][gs*gs];
    vector_matrix_row_multiply((double*)W_vector, (double*)J,
        (double*)tmp2, gs*gs, gs*gs*gs*gs);
    vector_matrix_row_multiply((double*)W_vector, (double*)eta,
        (double*)tmp3, gs*gs, gs*gs);
    cblas_dscal(gs*gs*gs*gs*gs*gs, 1.0/4/nx/ny*hx*hy, (double*)tmp2, 1);
    cblas_dscal(gs*gs*gs*gs, 1.0/4/nx/ny*hx*hy, (double*)tmp3, 1);
    // 定义M和B
    double* M = new double [gs*gs*gs*gs*nx*ny];
    double* B = new double [gs*gs*nx*ny];
    for(int i = 0; i < nx*ny; i++)
    {
        const double* k_local = k + gs*gs*i;
        const double* f_local = f + gs*gs*i;
        double tmpM[gs*gs][gs*gs*gs*gs];
        double tmpB[gs*gs][gs*gs];
        vector_matrix_row_multiply((double*)k_local, (double*)tmp2,
        (double*)tmpM, gs*gs, gs*gs*gs*gs); 
        vector_matrix_row_multiply((double*)f_local, (double*)tmp3,
        (double*)tmpB, gs*gs, gs*gs);
        double* sumM = new double[gs*gs*gs*gs];
        double* sumB = new double[gs*gs];
        matrix_column_sum((double*)tmpM, gs*gs, gs*gs*gs*gs, sumM);
        matrix_column_sum((double*)tmpB, gs*gs, gs*gs, sumB);
        for(int j = 0; j < gs*gs*gs*gs; j++)
        {
            M[i*gs*gs*gs*gs+j] = sumM[j];
        }
        for(int j = 0; j < gs*gs; j++)
        {
            B[i*gs*gs+j] = sumB[j];
        }
        delete [] sumM;
        delete [] sumB;
    }
    // 将s每行复制gs份拼凑到一起（张量积）
    int* c1 = new int [gs*gs*gs*gs*nx*ny];
    int* c2 = new int [gs*gs*gs*gs*nx*ny];
    for(int i = 0; i < nx*ny; i++)
    {
        for(int j = 0; j < gs*gs*gs*gs; j++)
        {
            c1[i*gs*gs*gs*gs+j] = s[i*gs*gs+j/(gs*gs)];
            c2[i*gs*gs*gs*gs+j] = s[i*gs*gs+j%(gs*gs)];
        }
    }
    double* b = new double [(nx*order+1)*(ny*order+1)];
    for(int i = 0; i < (nx*order+1)*(ny*order+1); i++)
    {
        b[i] = 0;
    }
    for(int i = 0; i < gs*gs*nx*ny; i++)
    {
        b[s[i]-1] += B[i];
    }
    delete [] s;
    
    
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
    delete [] c1;
    delete [] c2;
    delete [] M;
    delete [] c1_free;
    delete [] c2_free;
    delete [] M_free;
    delete [] B;
    delete [] b;
    delete [] pol_x;
    delete [] pol_y;
    delete [] t_x;
    delete [] t_y;
    delete [] W_vector;
    delete [] eta_x;
    delete [] eta_y;
    delete [] dex;
    delete [] dey;
    mkl_sparse_destroy(A);
    mkl_sparse_destroy(A_free);
    

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
    // print_vector("col_ind_free", col_ind_free, nnz_free);
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
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &sizefree, csr_val_free, 
        row_ptr_free, col_ind_free, nullptr, &nrhs, iparm, &msglvl, 
        nullptr, nullptr, &error);
    delete [] b_free;
    delete [] u_free;
    delete [] diri_val;
    delete [] diri_set;
    mkl_sparse_destroy(csr_A);
    mkl_sparse_destroy(csr_A_free);

    
}