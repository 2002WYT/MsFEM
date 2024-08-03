#include "head.h"
#include <iostream>
#include <iomanip>
using namespace std;
// 辅助函数用于初始化和打印矩阵
void print_vector(const char* name, double* A, int n) {
    cout << name << " =\n";
    for (int i = 0; i < n; i++) {
        cout << fixed << setprecision(4) << A[i] << " " << endl;
    }
}
void print_vector(const char* name, const double* A, int n) {
    cout << name << " =\n";
    for (int i = 0; i < n; i++) {
        cout << fixed << setprecision(4) << A[i] << " " << endl;
    }
}
// 辅助函数用于初始化和打印矩阵，进行函数重载，打印int型矩阵
void print_vector(const char* name, int* A, int n) {
    cout << name << " =\n";
    for (int i = 0; i < n; i++) {
        cout << fixed << setprecision(4) << A[i] << " " << endl;
    }
}
