#include "head.h"
void get_n_over(const int Nx, const int Ny, const int nx, const int ny, const int over, const int i_flag, 
    int& nx_over, int& ny_over)
{
    if (i_flag == 0) // 左下角
    {
        nx_over = nx + over;
        ny_over = ny + over;
    }
    else if (i_flag == Ny-1) // 左上角
    {
        nx_over = nx + over;
        ny_over = ny + over;
    }
    else if (i_flag == Ny*(Nx-1)) // 右下角
    {
        nx_over = nx + over;
        ny_over = ny + over;
    }
    else if (i_flag == Ny*Nx-1) // 右上角
    {
        nx_over = nx + over;
        ny_over = ny + over;
    }
    else if (i_flag > 0 && i_flag < Ny-1) // 左边
    {
        nx_over = nx + over;
        ny_over = ny + 2*over;
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == 0) // 下边
    {
        nx_over = nx + 2*over;
        ny_over = ny + over;
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == Ny-1) // 上边
    {
        nx_over = nx + 2*over;
        ny_over = ny + over;
    }
    else if (i_flag > Ny*(Nx-1) && i_flag < Ny*Nx-1) // 右边
    {
        nx_over = nx + over;
        ny_over = ny + 2*over;
    }
    else // 内部
    {
        nx_over = nx + 2*over;
        ny_over = ny + 2*over;
    }
    
}