// 将行向量乘以矩阵的每一行
#include "head.h"
#include <mkl.h>
void vector_matrix_row_multiply(const double* vec, const double* A, double* result, 
    const int m, const int n)
{
    mkl_domatcopy('R', 'N', m, n, 1.0, A, n, result, n);
    for(int i = 0; i < m; i++) {
        cblas_dscal(n, vec[i], result + i*n, 1);
    }
}