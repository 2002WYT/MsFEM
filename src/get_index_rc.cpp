#include "head.h"
// 用于获得振荡边界条件所需的单元边上的k和f的索引
void get_index_rc(const int Nx, const int Ny, const int nx, const int ny, const int nx_over, 
    const int ny_over, const int order, const int over, const int i_flag, int& index_c1,
    int& index_c2, int& index_r1, int& index_r2)
{
    const int gs = order + 1;
    if (i_flag == 0) // 左下角
    {
        index_c1 = 0;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = 0;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else if (i_flag == Ny-1) // 左上角
    {
        index_c1 = (Ny*ny-ny_over)*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = (Ny*ny-ny_over)*Nx*nx*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else if (i_flag == Ny*(Nx-1)) // 右下角
    {
        index_c1 = (Nx*nx-nx_over)*Ny*ny*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = (Nx*nx-nx_over)*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else if (i_flag == Ny*Nx-1) // 右上角
    {
        index_c1 = (Nx*nx-nx_over)*Ny*ny*gs + (Ny*ny-ny_over)*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = (Nx*nx-nx_over)*gs + (Ny*ny-ny_over)*Nx*nx*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else if (i_flag > 0 && i_flag < Ny-1) // 左边
    {
        index_c1 = (i_flag%Ny*ny-over)*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = (i_flag%Ny*ny-over)*Nx*nx*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == 0) // 下边
    {
        index_c1 = ((int)(i_flag/Ny)*nx-over)*Ny*ny*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = ((int)(i_flag/Ny)*nx-over)*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == Ny-1) // 上边
    {
        index_c1 = ((int)(i_flag/Ny)*nx-over)*Ny*ny*gs + (Ny*ny-ny_over)*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = ((int)(i_flag/Ny)*nx-over)*gs + (i_flag%Ny*ny-over)*Nx*nx*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else if (i_flag > Ny*(Nx-1) && i_flag < Nx*Ny-1) // 右边
    {
        index_c1 = (Nx*nx-nx_over)*Ny*ny*gs + (i_flag%Ny*ny-over)*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = (Nx*nx-nx_over)*gs + (i_flag%Ny*ny-over)*Nx*nx*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    else // 内部
    {
        index_c1 = ((int)(i_flag/Ny)*nx-over)*Ny*ny*gs + (i_flag%Ny*ny-over)*gs;
        index_c2 = index_c1 + nx_over*Ny*ny*gs;
        index_r1 = ((int)(i_flag/Ny)*nx-over)*gs + (i_flag%Ny*ny-over)*Nx*nx*gs;
        index_r2 = index_r1 + ny_over*Nx*nx*gs;
    }
    
    
    
    

}