#include "head.h"
void get_rightterm_free(int sizeA, double* b, int* diri_set, double* b_free)
{
    int cnt = 0;
    for(int i = 0; i < sizeA; i++)
    {
        if(diri_set[i] == -1)
        {
            continue;
        }
        else
        {
            b_free[cnt] = b[i];
            cnt++;
        }
    }
}