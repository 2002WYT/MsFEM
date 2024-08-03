#ifndef _HEAD_H
#define _HEAD_H
#include <mkl.h>
extern const double eps;
void gauss_interpolation(const int gs, double* inpol);
void gauss_coefficient(const int gs, double* coeff);
void print_matrix(const char* name, double* A, int m, int n);
void print_matrix(const char* name, int* A, int m, int n);
void matrix_to_vector(const double* A, double* A_vector, int m, int n);
void add_constant_to_vector(double* vec, const int size, const double constant);
void extract_column_as_vector(const double* A, const int N, const int M, const int k, double* vec);
double get_kappa(double x, double y, const int example);
double get_f(double x, double y, const int example);
double real_u(double x, double y, const int example);
// 得到高斯积分插值点坐标的函数
void get_total_coordinates_x(const int nx, const int ny, const int gs, const double x0, 
    const double hx, const double hy, 
    const double* pol_x, const double* pol_y, double* coorx);
void get_total_coordinates_y(const int nx, const int ny, const int gs, const double y0, 
    const double hx, const double hy, 
    const double* pol_x, const double* pol_y, double* coory);
// 向量每个元素与矩阵对应行乘积函数
void vector_matrix_row_multiply(const double* vec, const double* A, double* result, 
    const int m, const int n);
void matrix_row_sum(const double* A, const int m, const int n, double* result);
void matrix_column_sum(const double* A, const int m, const int n, double* result);
void print_vector(const char* name, double* A, int n);
void print_vector(const char* name, int* A, int n);
void print_vector(const char* name, const double* A, int n);
void mergeVectors(int* A, int* B, double* C, int& size);
void get_stiffmatrix_free(int nnz, int* c1, int* c2, double* M, int* c1_free, int* c2_free, double* M_free, 
    int* diri_set, int& nnz_free);
void get_rightterm_free(int sizeA, double* b, int* diri_set, double* b_free);
void get_meshgrid(const int nx, const int ny, const int order, const double x0, const double y0, 
    const double xl, const double yl, double* x, double* y);
void get_real_solution(double* meshx, double* meshy, double* ue, int sizemeshgrid, 
    const int example);
double get_relative_L2_norm(double* u, double* ue, int sizeu);
void solve_subfem2dho(const double x0, const double y0, const double xl, const double yl, 
    const int nx, const int ny, const int order, const double* bound_val, const double* inpol, 
    const double* coeff, const double* k, const double* f, double* u);
void get_stiffmatrix(const double x0, const double y0, const double xl, const double yl, 
    const int nx, const int ny, const int order, const double* inpol, 
    const double* coeff, const double* k, int* c1, int* c2, double* M);
void get_b(const double x0, const double y0, const double xl, const double yl, 
    const int nx, const int ny, const int order, const double* inpol, 
    const double* coeff, const double* f, double* b);
void calculate_Ab(const int nx, const int ny, const int order, const double* bound_val, int nnz, 
    int* c1, int* c2, double* M, double* b, double* u);
void get_block_coordinates_x(const int Nx, const int Ny, const int nx, const int ny, 
    const double x0, const int gs, 
    const double hx, const double hy, const double* pol_x, const double* pol_y, double* cocox);
void get_block_coordinates_y(const int Nx, const int Ny, const int nx, const int ny,
    const double y0, const int gs, 
    const double hx, const double hy, const double* pol_x, const double* pol_y, double* cocoy);
void get_x_column_first(const int Nx, const int Ny, const int nx, const int ny, const double x0, 
    const int gs, const double hx, const double* pol_x, double* x_column);
void get_x_row_first(const int Nx, const int Ny, const int nx, const int ny, const double x0, 
    const int gs, const double hx, const double* pol_x, double* x_row);
void get_y_column_first(const int Nx, const int Ny, const int nx, const int ny, const double y0, 
    const int gs, const double hy, const double* pol_y, double* y_column);
void get_y_row_first(const int Nx, const int Ny, const int nx, const int ny, const double y0, 
    const int gs, const double hy, const double* pol_y, double* y_row);
void get_s_fine(const int Nx, const int Ny, const int nx, const int ny, const int order, int* s_fine);
void get_poly_bound(const int order, const int nx, const int ny, double* pbd, const int flag);
void get_poly_term(const int order, const double* inpol, 
    const int nx, const int ny, const int flag, double* poly_term);
void sol1dho(const double x0, const double xl, const double* coeff, const double* inpol,
    const int n, const int order, const double* bound_val, const double* k, const double* f, 
    double* u);
void matrixcol_column_sum(const double* A, const int m, const int n, double* result);
void get_oscillatory_bound(const double x0, const double xl, const double y0, const double yl,
    const double* k_column1, const double* k_column2, const double* k_row1, const double* k_row2,
    const int nx, const int ny, const int order, const double* coeff, const double* inpol, 
    const int flag, double* osbd);
double get_relative_infinity_norm(double* u, double* ue, int sizeu);
void get_support_point(const int nx, const int ny, const int order, int* suppind);
void fix_basis_function(const int nx, const int ny, const int order, const double* R_local, 
    double* R_fix);
void print_matrix_col(const char* name, double* A, int m, int n);
void print_matrix_col(const char* name, int* A, int m, int n);
void get_over_domain(const double x0, const double x1, const double y0, const double y1, 
    const int Nx, const int Ny, const int nx, const int ny, const int i_flag, 
    const int over, double& left_over, double& right_over, double& bottom_over, double& top_over);
void get_n_over(const int Nx, const int Ny, const int nx, const int ny, const int over, 
    const int i_flag,  int& nx_over, int& ny_over);
void get_kf_domain(const int Nx, const int Ny, const int nx, const int ny, const int nx_over, 
    const int ny_over, const int order, const int over, const int i_flag, const double* k_total, 
    double* k_domain, const double* f_total, double* f_domain);
void get_index_rc(const int Nx, const int Ny, const int nx, const int ny, const int nx_over, 
    const int ny_over, const int order, const int over, const int i_flag, int& index_c1,
    int& index_c2, int& index_r1, int& index_r2);
void truncate_R(const int Nx, const int Ny, const int nx, const int ny, const int nx_over, 
    const int ny_over, const int order, const int over, const int i_flag, const double* R_over, 
    double* R_cut);
void get_oscillatory_bound_2(const double x0, const double xl, const double y0, const double yl,
                           const double *k_column1, const double *k_column2, const double *k_row1, const double *k_row2,
                           const int nx, const int ny, const int order, const double *coeff, const double *inpol,
                           double *osbd);
void get_poly_bound_2(const int order, const int nx, const int ny, double *pbd);
#endif
