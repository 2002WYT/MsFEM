#include "head.h"
#include <mkl.h>
void add_constant_to_vector(double* vec, const int size, const double constant) {
    // 创建一个与vec大小相同并且每个元素都等于constant的向量
    double* constant_vector = new double[size];
    for (int i = 0; i < size; ++i) {
        constant_vector[i] = constant;
    }
    for (int i = 0; i < size; ++i) {
        // 将vec中的每个元素加上constant
        vec[i] += constant_vector[i];
    }
    // 清理内存
    delete[] constant_vector;
}