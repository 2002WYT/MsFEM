#include "head.h"
void truncate_R(const int Nx, const int Ny, const int nx, const int ny, const int nx_over, 
    const int ny_over, const int order, const int over, const int i_flag, const double* R_over, 
    double* R_cut)
{
    const int gs = order + 1;
    if (i_flag == 0) // 左下角
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + i*(ny_over*order+1) + j];
                }
            }
        
        }
    }
    else if (i_flag == Ny-1) // 左上角
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + i*(ny_over*order+1) + j+over*order];
                }
            }
        
        }
    }
    else if (i_flag == Ny*(Nx-1)) // 右下角
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + (i+over*order)*(ny_over*order+1) + j];
                }
            }
        }
    }
    else if (i_flag == Ny*Nx-1) // 右上角
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + (i+over*order)*(ny_over*order+1) + j+over*order];
                }
            }
        }
    }
    else if (i_flag > 0 && i_flag < Ny-1) // 左边
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + i*(ny_over*order+1) + j+over*order];
                }
            }
        }
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == 0) // 下边
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + (i+over*order)*(ny_over*order+1) + j];
                }
            }
        }
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == Ny-1) // 上边
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + (i+over*order)*(ny_over*order+1) + j+over*order];
                }
            }
        }
    }
    else if (i_flag > Ny*(Nx-1) && i_flag < Nx*Ny-1) // 右边
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + (i+over*order)*(ny_over*order+1) + j+over*order];
                }
            }
        }
    }
    else // 内部
    {
        for (int q = 0; q < gs*gs; q++)
        {
            int qr_over = q*(nx_over*order+1)*(ny_over*order+1);
            int qr_cut = q*(nx*order+1)*(ny*order+1);
            for (int i = 0; i < nx*order+1; i++)
            {
                for (int j = 0; j < ny*order+1; j++)
                {
                    R_cut[qr_cut + i*(ny*order+1) + j] = 
                        R_over[qr_over + (i+over*order)*(ny_over*order+1) + j+over*order];
                }
            }
        }
    }
    
    
    
    

}