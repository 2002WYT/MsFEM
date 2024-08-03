#include "head.h"
void matrix_row_sum(const double* A, const int m, const int n, double* result)
{
    for(int i = 0; i < m; i++)
    {
        result[i] = 0;
        for(int j = 0; j < n; j++)
        {
            result[i] += A[i*n + j];
        }
    }
}