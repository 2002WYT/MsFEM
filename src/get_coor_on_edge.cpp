#include "head.h"
void get_x_column_first(const int Nx, const int Ny, const int nx, const int ny, const double x0, 
    const int gs, const double hx, const double* pol_x, double* x_column)
{
    for (int i = 0; i < Ny * ny * (Nx*nx+1) * gs; i++)
    {
        x_column[i] = (int)(i/(Ny*ny*gs)) * hx / Nx / nx;
        // 竖着计算在细网格边界的x坐标值
    }
}

void get_x_row_first(const int Nx, const int Ny, const int nx, const int ny, const double x0, 
    const int gs, const double hx, const double* pol_x, double* x_row)
{
    // 横着计算在细网格边界的x坐标值
    for (int i = 0; i < Ny * ny + 1; i++)
    {
        for (int j = 0; j < Nx * nx; j++)
        {
            double xj0 = x0 + j * hx / (Nx * nx);
            for (int k = 0; k < gs; k++)
            {
                x_row[i * (Nx * nx * gs) + j * gs + k] = pol_x[k] + xj0;
            }
            
        }
    }
}
void get_y_column_first(const int Nx, const int Ny, const int nx, const int ny, const double y0, 
    const int gs, const double hy, const double* pol_y, double* y_column)
{
    // 竖着计算在细网格边界的y坐标值
    for (int i = 0; i < Nx * nx + 1; i++)
    {
        for (int j = 0; j < Ny * ny; j++)
        {
            double yj0 = y0 + j * hy / (Ny * ny);
            for (int k = 0; k < gs; k++)
            {
                y_column[i * (Ny * ny * gs) + j * gs + k] = pol_y[k] + yj0;
            }
            
        }
    }
}
void get_y_row_first(const int Nx, const int Ny, const int nx, const int ny, const double y0, 
    const int gs, const double hy, const double* pol_y, double* y_row)
{
    for (int i = 0; i < Nx * nx * (Ny*ny+1) * gs; i++)
    {
        y_row[i] = (int)(i/(Nx*nx*gs)) * hy / Ny / ny;
        // 竖着计算在细网格边界的y坐标值
    }
}