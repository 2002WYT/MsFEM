#include "head.h"
void get_kf_domain(const int Nx, const int Ny, const int nx, const int ny, const int nx_over, 
    const int ny_over, const int order, const int over, const int i_flag, const double* k_total, 
    double* k_domain, const double* f_total, double* f_domain)
{
    const int gs = order + 1;
    if (i_flag == 0) // 左下角
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = i*Ny*ny*gs*gs + j*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else if (i_flag == Ny-1) // 左上角
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = i*Ny*ny*gs*gs + (Ny*ny-ny_over+j)*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else if (i_flag == Ny*(Nx-1)) // 右下角
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = (Nx*nx-nx_over+i)*Ny*ny*gs*gs + j*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else if (i_flag == Nx*Ny-1) // 右上角
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = (Nx*nx-nx_over+i)*Ny*ny*gs*gs + (Ny*ny-ny_over+j)*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else if (i_flag > 0 && i_flag < Ny-1) // 左边
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = i*Ny*ny*gs*gs + (i_flag*ny-over+j)*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == 0) // 下边
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = ((int)(i_flag/Ny)*nx-over+i)*Ny*ny*gs*gs + j*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else if (i_flag > Ny-1 && i_flag < Ny*(Nx-1) && i_flag%Ny == Ny-1) // 上边
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = ((int)(i_flag/Ny)*nx-over+i)*Ny*ny*gs*gs + 
                    (Ny*ny-ny_over+j)*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else if (i_flag > Ny*(Nx-1) && i_flag < Nx*Ny-1) // 右边
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = (Nx*nx-nx_over+i)*Ny*ny*gs*gs + (i_flag%Ny*ny-over+j)*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    else // 内部
    {
        for (int i = 0; i < nx_over; i++) // 先横向
        {
            for (int j = 0; j < ny_over; j++)
            {
                int index_domain = i*ny_over*gs*gs+j*gs*gs;
                int index_total = ((int)(i_flag/Ny)*nx-over+i)*Ny*ny*gs*gs + 
                    (i_flag%Ny*ny-over+j)*gs*gs;
                for (int q = 0; q < gs*gs; q++)
                {
                    k_domain[index_domain + q] = k_total[index_total + q];
                    f_domain[index_domain + q] = f_total[index_total + q];
                }
            }
        }
    }
    
}