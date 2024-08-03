// merge.cpp 把重复的坐标点和值合并
#include "head.h"
void mergeVectors(int* A, int* B, double* C, int& size) {
    int newSize = 0; // 新向量的长度
    bool* merged = new bool[size]; // 标记是否已经合并

    for (int i = 0; i < size; ++i) {
        if (!merged[i]) {
            for (int j = i + 1; j < size; ++j) {
                if (A[i] == A[j] && B[i] == B[j] && !merged[j]) {
                    C[i] += C[j];
                    merged[j] = true;
                }
            }
            A[newSize] = A[i];
            B[newSize] = B[i];
            C[newSize] = C[i];
            newSize++;
        }
    }
    delete [] merged; // 释放内存
    size = newSize; // 更新size为新的长度
}
// main.cpp
/* 示例矩阵
    int v1[6] = {1, 2, 1, 2, 1, 2};
    int v2[6] = {2, 3, 2, 3, 2, 5};
    double v3[6] = {9, 5, 8, 3, 7, 1};
    // 合并行
    int n = sizeof(v1) / sizeof(v1[0]);

    // 合并操作
    mergeVectors(v1, v2, v3, n);

    // 输出合并后的结果
    for (int i = 0; i < n; ++i) {
        cout << v1[i] << " " << v2[i] << " " << v3[i] << endl;
    }
*/