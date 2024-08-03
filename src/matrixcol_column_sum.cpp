#include "head.h"
void matrixcol_column_sum(const double* A, const int m, const int n, double* result)
{
    for(int j = 0; j < n; j++)
    {
        result[j] = 0;
        for(int i = 0; i < m; i++)
        {
            result[j] += A[j*m + i];
        }
    }
}