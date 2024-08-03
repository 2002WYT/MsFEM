#include "head.h"
void get_over_domain(const double x0, const double x1, const double y0, const double y1, 
    const int Nx, const int Ny, const int nx, const int ny, const int i_flag, 
    const int over, double& left_over, double& right_over, double& bottom_over, double& top_over)
{
    // Get the overation of the left, right, bottom, and top boundary
    // of the domain
    const double Lx = x1 - x0;
    const double Ly = y1 - y0;
    if (i_flag == 0) // 左下角
    {
        left_over = x0;
        right_over = left_over + Lx/Nx/nx*(nx+over);
        bottom_over = y0;
        top_over = bottom_over + Ly/Ny/ny*(ny+over);
    }
    else if (i_flag == Ny-1) // 左上角
    {
        left_over = x0;
        right_over = left_over + Lx/Nx/nx*(nx+over);
        bottom_over = y1 - Ly/Ny/ny*(ny+over);
        top_over = y1;
    }
    else if (i_flag == Ny*(Nx-1)) // 右下角
    {
        left_over = x1 - Lx/Nx/nx*(nx+over);
        right_over = x1;
        bottom_over = y0;
        top_over = bottom_over + Ly/Ny/ny*(ny+over);
    }
    else if (i_flag == Nx*Ny-1) // 右上角
    {
        left_over = x1 - Lx/Nx/nx*(nx+over);
        right_over = x1;
        bottom_over = y1 - Ly/Ny/ny*(ny+over);
        top_over = y1;
    }
    else if (i_flag > 0 && i_flag < Ny-1) // 左边
    {
        left_over = x0;
        right_over = left_over + Lx/Nx/nx*(nx+over);
        bottom_over = y0 + Ly/Ny/ny*(i_flag%Ny*ny-over);
        top_over = bottom_over + Ly/Ny/ny*(ny+over*2);
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == 0) // 下边
    {
        left_over = x0 + Lx/Nx/nx*((int)(i_flag/Ny)*nx-over);
        right_over = left_over + Lx/Nx/nx*(nx+over*2);
        bottom_over = y0;
        top_over = bottom_over + Ly/Ny/ny*(ny+over);
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == Ny-1) // 上边
    {
        left_over = x0 + Lx/Nx/nx*((int)(i_flag/Ny)*nx-over);
        right_over = left_over + Lx/Nx/nx*(nx+over*2);
        bottom_over = y1 - Ly/Ny/ny*(ny+over);
        top_over = y1;
    }
    else if (i_flag > Ny*(Nx-1) && i_flag < Nx*Ny-1) // 右边
    {
        left_over = x1 - Lx/Nx/nx*(nx+over);
        right_over = x1;
        bottom_over = y0 + Ly/Ny/ny*(i_flag%Ny*ny-over);
        top_over = bottom_over + Ly/Ny/ny*(ny+over*2);
    }
    else // 内部
    {
        left_over = x0 + Lx/Nx/nx*((int)(i_flag/Ny)*nx-over);
        right_over = left_over + Lx/Nx/nx*(nx+over*2);
        bottom_over = y0 + Ly/Ny/ny*(i_flag%Ny*ny-over);
        top_over = bottom_over + Ly/Ny/ny*(ny+over*2);
    }
}