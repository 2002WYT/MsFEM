// main function, including Petrov-Galerkin method and Oscillatory boundary condition
//      and oversampling method.
// type "ulimit -s unlimited" in ubuntu to avoid stack overflow!!!
// type "valgrind --leak-check=full ./main" to check memory leak
#include <iostream>
#include <mkl.h>
#include "head.h"
#include <mkl_spblas.h>
#include <cmath>
#include <chrono>
// #include <omp.h>
using namespace std;
int main()
{

    auto start_total = chrono::high_resolution_clock::now();

    double Time_osci = 0.0;
    double Time_poly = 0.0;
    cout << "This is MsFEM!" << endl;
    const int parallel_num = 16; // 并行线程数
    // bool ms = true; // 是否使用多尺度方法
    bool osci                    = 1;    // 是否使用振荡边界条件
    bool pg                      = 0;    // 是否使用Petrov-Galerkin方法
    bool real_solution_available = 0;    // 是否有真实解

    const int finemesh = 128;
    const int Nx       = 16;
    const int Ny       = 16;
    const int order    = 3;
    const int example  = 2;    // 选用例子
    if (finemesh % Nx != 0 || finemesh % Ny != 0)
    {
        cerr << "finemesh should be divisible by Nx and Ny" << endl;
        return -1;
    }
    const int nx   = finemesh / Nx;
    const int ny   = finemesh / Ny;
    const int over = ceil(1.0 / 4 * nx);    // 超采样单元数（细网格）
    double    x0 = 0.0, x1 = 1.0;
    double    y0 = 0.0, y1 = 1.0;
    if (order > 4)
    {
        cerr << "Order should be less than 5" << endl;
        return -1;
    }
    const int gs    = order + 1;
    double   *inpol = (double *)mkl_malloc(gs * sizeof(double), 64);
    for (int i = 0; i < gs; i++)
    {
        inpol[i] = 0.0;
    }
    gauss_interpolation(gs, inpol);
    double *coeff = (double *)mkl_malloc(gs * sizeof(double), 64);
    for (int i = 0; i < gs; i++)
    {
        coeff[i] = 0.0;
    }
    gauss_coefficient(gs, coeff);
    double *bound_val = (double *)mkl_malloc((Nx * order + 1) * (Ny * order + 1) * sizeof(double), 64);
    for (int i = 0; i < (Nx * order + 1) * (Ny * order + 1); i++)
    {
        bound_val[i] = 0;
    }

    double    Lx       = x1 - x0;
    double    Ly       = y1 - y0;
    const int sizeA    = (Nx * order + 1) * (Ny * order + 1);
    const int sizefine = (finemesh * order + 1) * (finemesh * order + 1);
    double   *pol_x    = (double *)mkl_malloc(gs * sizeof(double), 64);
    double   *pol_y    = (double *)mkl_malloc(gs * sizeof(double), 64);
    for (int i = 0; i < gs; i++)
    {
        pol_x[i] = inpol[i];
        pol_y[i] = inpol[i];
    }
    // 利用高斯插值点计算网格内的所需高斯积分插值点
    add_constant_to_vector(pol_x, gs, 1.0);
    add_constant_to_vector(pol_y, gs, 1.0);
    cblas_dscal(gs, 1.0 / 2 * Lx / finemesh, (double *)pol_x, 1);
    cblas_dscal(gs, 1.0 / 2 * Ly / finemesh, (double *)pol_y, 1);
    double *cocox = (double *)mkl_malloc(Nx * Ny * nx * ny * gs * gs * sizeof(double),
                                         64);
    double *cocoy = (double *)mkl_malloc(Nx * Ny * nx * ny * gs * gs * sizeof(double),
                                         64);
    // 按照粗网格顺序排列
    get_block_coordinates_x(Nx, Ny, nx, ny, x0, gs, Lx, Ly, pol_x, pol_y, cocox);
    get_block_coordinates_y(Nx, Ny, nx, ny, y0, gs, Lx, Ly, pol_x, pol_y, cocoy);
    double *k_block = (double *)mkl_malloc(Nx * Ny * nx * ny * gs * gs * sizeof(double),
                                           64);
    double *f_block = (double *)mkl_malloc(Nx * Ny * nx * ny * gs * gs * sizeof(double),
                                           64);
    for (int i = 0; i < Nx * Ny * nx * ny * gs * gs; i++)
    {
        k_block[i] = get_kappa(cocox[i], cocoy[i], example);
        f_block[i] = get_f(cocox[i], cocoy[i], example);
    }
    double *x_column = (double *)mkl_malloc(Ny * ny * gs * (Nx * nx + 1) * sizeof(double),
                                            64);
    double *y_column = (double *)mkl_malloc(Ny * ny * gs * (Nx * nx + 1) * sizeof(double),
                                            64);
    double *x_row    = (double *)mkl_malloc(Nx * nx * gs * (Ny * ny + 1) * sizeof(double),
                                            64);
    double *y_row    = (double *)mkl_malloc(Nx * nx * gs * (Ny * ny + 1) * sizeof(double),
                                            64);
    get_x_column_first(Nx, Ny, nx, ny, x0, gs, Lx, pol_x, x_column);
    get_y_column_first(Nx, Ny, nx, ny, y0, gs, Ly, pol_y, y_column);
    get_x_row_first(Nx, Ny, nx, ny, x0, gs, Lx, pol_x, x_row);
    get_y_row_first(Nx, Ny, nx, ny, y0, gs, Ly, pol_y, y_row);
    double *k_column = (double *)mkl_malloc(Ny * ny * gs * (Nx * nx + 1) * sizeof(double),
                                            64);
    double *f_column = (double *)mkl_malloc(Ny * ny * gs * (Nx * nx + 1) * sizeof(double),
                                            64);
    double *k_row    = (double *)mkl_malloc(Nx * nx * gs * (Ny * ny + 1) * sizeof(double),
                                            64);
    double *f_row    = (double *)mkl_malloc(Nx * nx * gs * (Ny * ny + 1) * sizeof(double),
                                            64);
    for (int i = 0; i < Ny * ny * gs * (Nx * nx + 1); i++)
    {
        k_column[i] = get_kappa(x_column[i], y_column[i], example);
        f_column[i] = get_f(x_column[i], y_column[i], example);
    }
    for (int i = 0; i < Nx * nx * gs * (Ny * ny + 1); i++)
    {
        k_row[i] = get_kappa(x_row[i], y_row[i], example);
        f_row[i] = get_f(x_row[i], y_row[i], example);
    }
    double *x_total = (double *)mkl_malloc(gs * gs * Nx * Ny * nx * ny * sizeof(double),
                                           64);
    get_total_coordinates_x(Nx * nx, Ny * ny, gs, x0, Lx, Ly, pol_x, pol_y,
                            x_total);
    double *y_total = (double *)mkl_malloc(gs * gs * Nx * Ny * nx * ny * sizeof(double),
                                           64);
    get_total_coordinates_y(Nx * nx, Ny * ny, gs, y0, Lx, Ly, pol_x, pol_y,
                            y_total);
    double *k_total = (double *)mkl_malloc(gs * gs * nx * ny * Nx * Ny * sizeof(double),
                                           64);
    double *f_total = (double *)mkl_malloc(gs * gs * nx * ny * Nx * Ny * sizeof(double),
                                           64);
    for (int i = 0; i < gs * gs * nx * ny * Nx * Ny; i++)
    {
        k_total[i] = get_kappa(x_total[i], y_total[i], example);
        f_total[i] = get_f(x_total[i], y_total[i], example);
    }    // 按照一列一列记录细网格
    int *s      = (int *)mkl_malloc(Nx * Ny * gs * gs * sizeof(int), 64);
    int *s_fine = (int *)mkl_malloc(Nx * Ny * (nx * order + 1) * (ny * order + 1) * sizeof(int), 64);
    get_s_fine(Nx, Ny, 1, 1, order, s);           // 粗网格格点标号，用作组装A和b
    get_s_fine(Nx, Ny, nx, ny, order, s_fine);    // 细网格格点标号，用作最后投影至细网格解
    double *t_x = (double *)mkl_malloc(gs * sizeof(double), 64);
    double *t_y = (double *)mkl_malloc(gs * sizeof(double), 64);
    for (int i = 0; i < gs; i++)
    {
        t_x[i] = Lx / Nx / order * i;
        t_y[i] = Ly / Ny / order * i;
    }
    double *projx = (double *)mkl_malloc((nx * order + 1) * gs * sizeof(double), 64);
    double *projy = (double *)mkl_malloc((ny * order + 1) * gs * sizeof(double), 64);
    for (int i = 0; i < nx * order + 1; i++)
    {
        projx[i] = Lx / Nx / nx / order * i;
    }
    for (int i = 0; i < ny * order + 1; i++)
    {
        projy[i] = Ly / Ny / ny / order * i;
    }
    double *eta_x2 = (double *)mkl_malloc((nx * order + 1) * gs * sizeof(double), 64);
    double *eta_y2 = (double *)mkl_malloc((ny * order + 1) * gs * sizeof(double), 64);
    for (int i = 0; i < (nx * order + 1) * gs; i++)
    {
        eta_x2[i] = 1.0;
    }
    for (int i = 0; i < (ny * order + 1) * gs; i++)
    {
        eta_y2[i] = 1.0;
    }    // 列优先
    for (int i = 0; i < gs; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            if (j != i)
            {
                for (int k = 0; k < nx * order + 1; k++)
                {
                    eta_x2[i * (nx * order + 1) + k] *= (projx[k] - t_x[j]) / (t_x[i] - t_x[j]);
                }
                for (int k = 0; k < ny * order + 1; k++)
                {
                    eta_y2[i * (ny * order + 1) + k] *= (projy[k] - t_y[j]) / (t_y[i] - t_y[j]);
                }
            }
        }
    }
    double *R_proj = (double *)mkl_malloc((ny * order + 1) * (nx * order + 1) * gs * gs * sizeof(double), 64);
    // 列优先的多项式基函数矩阵
    for (int i = 0; i < gs; i++)
    {
        for (int j = 0; j < gs; j++)
        {
            int ijr = i * gs * (ny * order + 1) * (nx * order + 1) + j * (ny * order + 1) * (nx * order + 1);
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                        ny * order + 1, nx * order + 1, 1, 1.0, eta_y2 + j * (ny * order + 1),
                        ny * order + 1, eta_x2 + i * (nx * order + 1), nx * order + 1, 0.0,
                        R_proj + ijr, ny * order + 1);
        }
    }


    int *flag_index = (int *)mkl_malloc(gs * gs * sizeof(int), 64);
    for (int i = 0; i < gs * gs; i++)
    {
        if (i % gs != order && i % gs != 0 && i > gs - 1 && i < gs * gs - gs)
            flag_index[i] = 0;
        else
            flag_index[i] = 1;    // 支撑点在单元内部的函数序号标0，否则标1
    }

    double *R = (double *)mkl_malloc((ny * order + 1) * (nx * order + 1) * gs * gs * Ny * Nx * sizeof(double),
                                     64);
    // 存储每个单元的基函数
    double *M = (double *)mkl_malloc(gs * gs * gs * gs * Nx * Ny * sizeof(double), 64);
    double *B = (double *)mkl_malloc(gs * gs * Nx * Ny * sizeof(double), 64);

    // 接下来开始遍历各个粗网格单元计算粗尺度基函数。
    #pragma omp parallel for num_threads(parallel_num)
    for (int i = 0; i < Nx * Ny; i++)
    {
        double left_loc   = floor(i / Ny) * Lx / Nx;
        double right_loc  = left_loc + Lx / Nx;
        double bottom_loc = (i % Ny) * Ly / Ny;
        double top_loc    = bottom_loc + Ly / Ny;
        double left_over, right_over, bottom_over, top_over;
        get_over_domain(x0, x1, y0, y1, Nx, Ny, nx, ny, i, over,
                        left_over, right_over, bottom_over, top_over);
        int nx_over, ny_over;
        get_n_over(Nx, Ny, nx, ny, over, i, nx_over, ny_over);
        double *k_domain = (double *)mkl_malloc(gs * gs * nx_over * ny_over * sizeof(double),
                                                64);
        double *f_domain = (double *)mkl_malloc(gs * gs * nx_over * ny_over * sizeof(double),
                                                64);
        get_kf_domain(Nx, Ny, nx, ny, nx_over, ny_over, order, over, i, k_total,
                      k_domain, f_total, f_domain);
        int indexc1, indexc2, indexr1, indexr2;
        get_index_rc(Nx, Ny, nx, ny, nx_over, ny_over, order, over, i, indexc1,
                     indexc2, indexr1, indexr2);
        const int sizebd_over = 2 * nx_over * order + 2 * ny_over * order;
        int       sizeAsmall  = (nx * order + 1) * (ny * order + 1);
        const int nnzsmall    = gs * gs * gs * gs * nx * ny;
        int      *c1small     = (int *)mkl_malloc(nnzsmall * sizeof(int), 64);
        int      *c2small     = (int *)mkl_malloc(nnzsmall * sizeof(int), 64);
        double   *Msmall      = (double *)mkl_malloc(nnzsmall * sizeof(double), 64);
        double   *b_small     = (double *)mkl_malloc(sizeAsmall * sizeof(double), 64);
        // A_small用作日后R_test^T * A_small * R_test，大小正常
        get_stiffmatrix(left_loc, bottom_loc, right_loc, top_loc, nx, ny,
                        order, inpol,
                        coeff, k_block + i * nx * ny * gs * gs, c1small, c2small, Msmall);
        sparse_matrix_t A_small;
        mkl_sparse_d_create_coo(&A_small, SPARSE_INDEX_BASE_ONE, sizeAsmall,
                                sizeAsmall, nnzsmall, c1small, c2small, Msmall);
        sparse_matrix_t csr_A_small;
        mkl_sparse_convert_csr(A_small, SPARSE_OPERATION_NON_TRANSPOSE,
                               &csr_A_small);

        int       sizeAover = (nx_over * order + 1) * (ny_over * order + 1);
        const int nnzover   = gs * gs * gs * gs * nx_over * ny_over;
        int      *c1over    = (int *)mkl_malloc(nnzover * sizeof(int), 64);
        int      *c2over    = (int *)mkl_malloc(nnzover * sizeof(int), 64);
        double   *Mover     = (double *)mkl_malloc(nnzover * sizeof(double), 64);
        // A_over用作求超采样范围内基函数，大小为超采样区域
        get_stiffmatrix(left_over, bottom_over, right_over, top_over, nx_over,
                        ny_over, order, inpol, coeff, k_domain, c1over, c2over, Mover);

        double *R_over = (double *)mkl_malloc(sizeAover * gs * gs * sizeof(double), 64);
        double *R_cut  = (double *)mkl_malloc(sizeAsmall * gs * gs * sizeof(double), 64);
        if (osci == 1)    // 振荡边界条件
        {
            double *osbd = (double *)mkl_malloc(sizebd_over * gs * gs * sizeof(double), 64);

            get_oscillatory_bound_2(left_over, right_over, bottom_over, top_over,
                                    k_column + indexc1, k_column + indexc2, k_row + indexr1, k_row + indexr2,
                                    nx_over, ny_over, order, coeff, inpol, osbd);

            double *zeroterm = (double *)mkl_malloc(sizeAover * sizeof(double), 64);
            double *zerobd   = (double *)mkl_malloc(sizebd_over * sizeof(double), 64);
            for (int j = 0; j < gs * gs; j++)
            {
                if (flag_index[j] == 1)    // 基函数支撑点在边界
                {
                    for (int q = 0; q < sizeAover; q++)
                    {
                        zeroterm[q] = 0.0;
                    }
                    int  jr         = j * sizeAover;
                    
                    calculate_Ab(nx_over, ny_over, order, osbd + j * sizebd_over, nnzover,
                                 c1over, c2over, Mover, zeroterm, R_over + jr);
                    
                }
                else    // 基函数支撑点在内部
                {
                    for (int q = 0; q < sizebd_over; q++)
                    {
                        zerobd[q] = 0.0;
                    }
                    auto start_osci = chrono::high_resolution_clock::now();
                    double *poly_term = (double *)mkl_malloc(gs * gs * nx_over * ny_over
                                                                 * sizeof(double),
                                                             64);
                    for (int q = 0; q < gs * gs * nx_over * ny_over; q++)
                    {
                        poly_term[q] = 1.0;
                    }
                    get_poly_term(order, inpol, nx_over, ny_over, j, poly_term);
                    // 得到粗网格内的多项式基函数用作以后细尺度有限元的右端项。
                    auto end_osci      = chrono::high_resolution_clock::now();
                    auto duration_poly = chrono::duration_cast< chrono::microseconds >(end_osci - start_osci);
                    Time_osci += (double)duration_poly.count() / 1e6;
                    double *b_small2 = (double *)mkl_malloc(sizeAover * sizeof(double), 64);
                    
                    get_b(left_over, bottom_over, right_over, top_over, nx_over,
                          ny_over, order, inpol,
                          coeff, poly_term, b_small2);
                          
                    calculate_Ab(nx_over, ny_over, order, zerobd, nnzover,
                                 c1over, c2over, Mover, b_small2, R_over + j * sizeAover);
                    
                    mkl_free(b_small2);
                    mkl_free(poly_term);
                }
            }
            mkl_free(zeroterm);
            mkl_free(zerobd);
            mkl_free(osbd);
        }
        else if (osci == 0)   // 多项式边界条件
        {
            auto    start_poly = chrono::high_resolution_clock::now();
            double *pybd       = (double *)mkl_malloc(sizebd_over * gs * gs * sizeof(double), 64);
            get_poly_bound_2(order, nx_over, ny_over, pybd);
            double *zeroterm = (double *)mkl_malloc(sizeAover * sizeof(double), 64);
            double *zerobd   = (double *)mkl_malloc(sizebd_over * sizeof(double), 64);
            for (int j = 0; j < gs * gs; j++)
            {
                if (flag_index[j] == 1)
                {
                    for (int q = 0; q < sizeAover; q++)
                    {
                        zeroterm[q] = 0.0;
                    }
                    int jr = j * sizeAover;
                    calculate_Ab(nx_over, ny_over, order, pybd + j * sizebd_over, nnzover,
                                 c1over, c2over, Mover, zeroterm, R_over + jr);
                }
                else if (flag_index[j] == 0)
                {
                    for (int q = 0; q < sizebd_over; q++)
                    {
                        zerobd[q] = 0.0;
                    }
                    double *poly_term = (double *)mkl_malloc(gs * gs * nx_over * ny_over
                                                                 * sizeof(double),
                                                             64);
                    for (int q = 0; q < gs * gs * nx_over * ny_over; q++)
                    {
                        poly_term[q] = 1.0;
                    }
                    get_poly_term(order, inpol, nx_over, ny_over, j, poly_term);
                    // 得到粗网格内的多项式基函数用作以后细尺度有限元的右端项。
                    double *b_small2 = (double *)mkl_malloc(sizeAover * sizeof(double), 64);
                    get_b(left_over, bottom_over, right_over, top_over, nx_over,
                          ny_over, order, inpol,
                          coeff, poly_term, b_small2);
                    calculate_Ab(nx_over, ny_over, order, zerobd, nnzover,
                                 c1over, c2over, Mover, b_small2, R_over + j * sizeAover);
                    mkl_free(b_small2);
                    mkl_free(poly_term);
                }
            }
            mkl_free(zeroterm);
            mkl_free(zerobd);
            mkl_free(pybd);
            auto end_poly      = chrono::high_resolution_clock::now();
            auto duration_poly = chrono::duration_cast< chrono::microseconds >(end_poly - start_poly);
            Time_poly += (double)duration_poly.count() / 1e6;
        }
        truncate_R(Nx, Ny, nx, ny, nx_over, ny_over, order, over, i, R_over, R_cut);

        double *R_fix = (double *)mkl_malloc(gs * gs * sizeAsmall * sizeof(double), 64);

        fix_basis_function(nx, ny, order, R_cut, R_fix);
        // 将基函数线性组合使之在对应点处值为1和0


        // 将R_fix记录到R中
        for (int j = 0; j < gs * gs * sizeAsmall; j++)
        {
            R[i * gs * gs * sizeAsmall + j] = R_fix[j];
        }

        get_b(left_loc, bottom_loc, right_loc, top_loc, nx, ny, order, inpol,
              coeff, f_block + i * nx * ny * gs * gs, b_small);
        double *R_test = (double *)mkl_malloc(gs * gs * sizeAsmall * sizeof(double), 64);
        if (pg == true)    // Petrov-Galerkin方法
        {
            for (int q = 0; q < gs * gs * sizeAsmall; q++)
            {
                R_test[q] = R_proj[q];
            }
            // print_matrix_col("R_test", R_test, sizeAsmall, gs*gs);
        }
        else    // Galerkin方法
        {
            for (int q = 0; q < gs * gs * sizeAsmall; q++)
            {
                R_test[q] = R_fix[q];
            }
        }
        double *R_trial = (double *)mkl_malloc(gs * gs * sizeAsmall * sizeof(double), 64);
        for (int q = 0; q < gs * gs * sizeAsmall; q++)
        {
            R_trial[q] = R_fix[q];
        }

        double *tmpRA = (double *)mkl_malloc(gs * gs * sizeAsmall * sizeof(double), 64);
        for (int q = 0; q < gs * gs * sizeAsmall; q++)
        {
            tmpRA[q] = 0.0;
        }
        matrix_descr descr1;
        descr1.type = SPARSE_MATRIX_TYPE_GENERAL;
        mkl_sparse_d_mm(SPARSE_OPERATION_TRANSPOSE, 1.0, csr_A_small,
                        descr1, SPARSE_LAYOUT_COLUMN_MAJOR, R_test, gs * gs,
                        sizeAsmall, 0.0, tmpRA, sizeAsmall);
        // tmpRA = A_small^T * R_test
        // 2024.07.23 WYT由于用错了columns: gs*gs, 导致内存一直泄露找不出问题，现已修正
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, gs * gs,
                    gs * gs, sizeAsmall, 1.0, tmpRA, sizeAsmall,
                    R_trial, sizeAsmall, 0.0, M + i * gs * gs * gs * gs, gs * gs);
        // R_test^T * A_small * R_trial

        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, gs * gs,
                    1, sizeAsmall, 1.0, R_test, sizeAsmall, b_small,
                    sizeAsmall, 0.0, B + i * gs * gs, gs * gs);
        // R_test^T * b_small

        mkl_free(k_domain);
        mkl_free(f_domain);
        mkl_free(R_cut);
        mkl_free(R_test);
        mkl_free(R_trial);
        mkl_free(c1small);
        mkl_free(c2small);
        mkl_free(Msmall);
        mkl_free(b_small);
        mkl_free(tmpRA);
        mkl_free(R_fix);
        mkl_sparse_destroy(A_small);
        mkl_sparse_destroy(csr_A_small);
        mkl_free(c1over);
        mkl_free(c2over);
        mkl_free(Mover);
        mkl_free(R_over);
    }

    int *c1 = (int *)mkl_malloc(gs * gs * gs * gs * Nx * Ny * sizeof(int), 64);
    int *c2 = (int *)mkl_malloc(gs * gs * gs * gs * Nx * Ny * sizeof(int), 64);
    for (int i = 0; i < Nx * Ny; i++)
    {
        for (int j = 0; j < gs * gs * gs * gs; j++)
        {
            c1[i * gs * gs * gs * gs + j] = s[i * gs * gs + j / (gs * gs)];
            c2[i * gs * gs * gs * gs + j] = s[i * gs * gs + j % (gs * gs)];
        }
    }
    double *b = (double *)mkl_malloc((Nx * order + 1) * (Ny * order + 1) * sizeof(double),
                                     64);
    for (int i = 0; i < (Nx * order + 1) * (Ny * order + 1); i++)
    {
        b[i] = 0;
    }
    for (int i = 0; i < gs * gs * Nx * Ny; i++)
    {
        b[s[i] - 1] += B[i];
    }
    double *uc = (double *)mkl_malloc(sizeA * sizeof(double), 64);
    // print_matrix_col("M", M, gs*gs*gs*gs, Nx*Ny);

    calculate_Ab(Nx, Ny, order, bound_val, gs * gs * gs * gs * Nx * Ny, c1, c2, M, b,
                 uc);

    // 以下计算FEM粗尺度解
    double *pol_xc = (double *)mkl_malloc(gs * sizeof(double), 64);
    double *pol_yc = (double *)mkl_malloc(gs * sizeof(double), 64);
    for (int i = 0; i < gs; i++)
    {
        pol_xc[i] = inpol[i];
        pol_yc[i] = inpol[i];
    }
    // 利用高斯插值点计算网格内的所需高斯积分插值点
    add_constant_to_vector(pol_xc, gs, 1.0);
    add_constant_to_vector(pol_yc, gs, 1.0);
    cblas_dscal(gs, 1.0 / 2 * Lx / Nx, pol_xc, 1);
    cblas_dscal(gs, 1.0 / 2 * Ly / Ny, pol_yc, 1);
    double *xc = (double *)mkl_malloc(gs * gs * Nx * Ny * sizeof(double), 64);
    get_total_coordinates_x(Nx, Ny, gs, x0, Lx, Ly, pol_xc, pol_yc,
                            xc);
    double *yc = (double *)mkl_malloc(gs * gs * Nx * Ny * sizeof(double), 64);
    get_total_coordinates_y(Nx, Ny, gs, y0, Lx, Ly, pol_xc, pol_yc,
                            yc);
    double *kc = (double *)mkl_malloc(gs * gs * Nx * Ny * sizeof(double), 64);
    double *fc = (double *)mkl_malloc(gs * gs * Nx * Ny * sizeof(double), 64);
    for (int i = 0; i < gs * gs * Nx * Ny; i++)
    {
        kc[i] = get_kappa(xc[i], yc[i], example);
        fc[i] = get_f(xc[i], yc[i], example);
    }
    double *u2c = (double *)mkl_malloc(sizeA * sizeof(double), 64);
    solve_subfem2dho(x0, y0, x1, y1, Nx, Ny, order, bound_val, inpol, coeff,
                     kc, fc, u2c);    // FEM粗尺度解
    // 打印粗尺度解
    // print_vector("uc", uc, sizeA);
    if (real_solution_available)
    {
        double *mesLx = (double *)mkl_malloc(sizeA * sizeof(double), 64);
        double *mesLy = (double *)mkl_malloc(sizeA * sizeof(double), 64);
        get_meshgrid(Nx, Ny, order, x0, y0, x1, y1, mesLx, mesLy);
        double *uce = (double *)mkl_malloc(sizeA * sizeof(double), 64);
        get_real_solution(mesLx, mesLy, uce, sizeA, example);
        double L2_normuc = get_relative_L2_norm(uc, uce, sizeA);
        cout << "MSFEM粗尺度解L2 norm: " << L2_normuc << endl;
        double L2_normu2c = get_relative_L2_norm(u2c, uce, sizeA);
        cout << "FEM粗尺度解L2 norm: " << L2_normu2c << endl;
        mkl_free(mesLx);
        mkl_free(mesLy);
        mkl_free(uce);
    }
    // MSFEM细尺度解
    double *uf = (double *)mkl_malloc(sizefine * sizeof(double), 64);
    // FEM细尺度解
    double *u2f = (double *)mkl_malloc(sizefine * sizeof(double), 64);
    for (int i = 0; i < Nx * Ny; i++)
    {
        double *T = (double *)mkl_malloc(gs * gs * sizeof(double), 64);
        // 存储一个粗网格内支撑点的解的值
        double *T2    = (double *)mkl_malloc(gs * gs * sizeof(double), 64);
        double *tmpf  = (double *)mkl_malloc((nx * order + 1) * (ny * order + 1) * sizeof(double), 64);
        double *tmpf2 = (double *)mkl_malloc((nx * order + 1) * (ny * order + 1) * sizeof(double), 64);
        for (int j = 0; j < gs * gs; j++)
        {
            T[j]  = uc[s[i * gs * gs + j] - 1];
            T2[j] = u2c[s[i * gs * gs + j] - 1];
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (nx * order + 1) * (ny * order + 1),
                    1, gs * gs, 1.0, R + i * gs * gs * (nx * order + 1) * (ny * order + 1),
                    (nx * order + 1) * (ny * order + 1), T, gs * gs, 0.0, tmpf,
                    (nx * order + 1) * (ny * order + 1));
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (nx * order + 1) * (ny * order + 1),
                    1, gs * gs, 1.0, R_proj,
                    (nx * order + 1) * (ny * order + 1), T2, gs * gs, 0.0, tmpf2,
                    (nx * order + 1) * (ny * order + 1));
        for (int j = 0; j < (nx * order + 1) * (ny * order + 1); j++)
        {
            uf[s_fine[i * (nx * order + 1) * (ny * order + 1) + j] - 1]  = tmpf[j];
            u2f[s_fine[i * (nx * order + 1) * (ny * order + 1) + j] - 1] = tmpf2[j];
        }


        mkl_free(T);
        mkl_free(T2);
        mkl_free(tmpf);
        mkl_free(tmpf2);
    }
    if (real_solution_available)
    {
        double *meshxf = (double *)mkl_malloc(sizefine * sizeof(double), 64);
        double *meshyf = (double *)mkl_malloc(sizefine * sizeof(double), 64);
        get_meshgrid(Nx * nx, Ny * ny, order, x0, y0, x1, y1, meshxf, meshyf);
        double *ufe = (double *)mkl_malloc(sizefine * sizeof(double), 64);
        get_real_solution(meshxf, meshyf, ufe, sizefine, example);
        double L2_normuf = get_relative_L2_norm(uf, ufe, sizefine);
        cout << "MSFEM细尺度解L2 norm: " << L2_normuf << endl;
        double L2_normu2f = get_relative_L2_norm(u2f, ufe, sizefine);
        cout << "FEM细尺度解L2 norm: " << L2_normu2f << endl;
        mkl_free(meshxf);
        mkl_free(meshyf);
        mkl_free(ufe);
    }
    else
    {
        double *bound_valf = (double *)mkl_malloc((Nx * nx * order + 1) * (Ny * ny * order + 1) * sizeof(double), 64);
        for (int i = 0; i < (Nx * nx * order + 1) * (Ny * ny * order + 1); i++)
        {
            bound_valf[i] = 0;
        }
        double *uref = (double *)mkl_malloc(sizefine * sizeof(double), 64);
        solve_subfem2dho(x0, y0, x1, y1, Nx * nx, Ny * ny, order, bound_valf,
                         inpol, coeff, k_total, f_total, uref);
        double L2_normuf = get_relative_L2_norm(uf, uref, sizefine);
        cout << "MSFEM细尺度解L2 norm: " << L2_normuf << endl;
        double L2_normu2f = get_relative_L2_norm(u2f, uref, sizefine);
        cout << "FEM细尺度解L2 norm: " << L2_normu2f << endl;
        mkl_free(bound_valf);
        mkl_free(uref);
    }


    // 释放内存
    mkl_free(inpol);
    mkl_free(coeff);
    mkl_free(bound_val);
    mkl_free(pol_x);
    mkl_free(pol_y);
    mkl_free(cocox);
    mkl_free(cocoy);
    mkl_free(k_block);
    mkl_free(f_block);
    mkl_free(x_column);
    mkl_free(y_column);
    mkl_free(k_column);
    mkl_free(f_column);
    mkl_free(x_row);
    mkl_free(y_row);
    mkl_free(k_row);
    mkl_free(f_row);
    mkl_free(x_total);
    mkl_free(y_total);
    mkl_free(k_total);
    mkl_free(f_total);
    mkl_free(s);
    mkl_free(s_fine);
    mkl_free(t_x);
    mkl_free(t_y);
    mkl_free(projx);
    mkl_free(projy);
    mkl_free(R_proj);
    mkl_free(flag_index);
    mkl_free(R);
    mkl_free(M);
    mkl_free(B);
    mkl_free(eta_x2);
    mkl_free(eta_y2);
    mkl_free(c1);
    mkl_free(c2);
    mkl_free(b);
    mkl_free(uc);
    mkl_free(pol_xc);
    mkl_free(pol_yc);
    mkl_free(xc);
    mkl_free(yc);
    mkl_free(kc);
    mkl_free(fc);
    mkl_free(u2c);
    mkl_free(uf);
    mkl_free(u2f);


    auto end_total      = chrono::high_resolution_clock::now();
    auto duration_total = chrono::duration_cast< chrono::microseconds >(end_total - start_total);

    cout << "Time taken by osci: "
         << Time_osci << " seconds" << endl;
    cout << "Total time taken: "
         << (double)duration_total.count() / 1e6 << " seconds" << endl;

    return 0;
}