#include "head.h"
#include <iostream>
#include <iomanip>
using namespace std;
// 辅助函数用于初始化和打印矩阵
void print_matrix(const char* name, double* A, int m, int n) {
    int width = 10;
    cout << name << " =\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << setw(width) << A[i * n + j] << " ";
        }
        cout << endl;
    }
}
// 辅助函数用于初始化和打印矩阵，进行函数重载，打印int型矩阵
void print_matrix(const char* name, int* A, int m, int n) {
    int width = 10;
    cout << name << " =\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << setw(width) << A[i * n + j] << " ";
        }
        cout << endl;
    }
}
