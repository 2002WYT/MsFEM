#include "head.h"
void extract_column_as_vector(const double* A, const int N, const int M, const int k, double* vec) {
    // A是N*M矩阵，vec是长度为N的向量
    // 提取第k列到vec中
    for (int i = 0; i < N; ++i) {
        vec[i] = A[i * M + k];
    }
}