#include "head.h"
// 把矩阵变为行向量，按照行优先存储。
void matrix_to_vector(const double* A, double* A_vector, int m, int n) {
    for (int col = 0; col < n; col++) {
        for (int row = 0; row < m; row++) {
            A_vector[col * m + row] = A[row * n + col];
        }
    }
}