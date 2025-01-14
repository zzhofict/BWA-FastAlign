/*********************************************************************************************
    Description: data and files for debugging
    Copyright : All right reserved by .

    Author : 
    Date : 2024/04/09
***********************************************************************************************/
#include <stdio.h>
#include "debug.h"

FILE *gf[4] = {0},
     *gfq[4] = {0},
     *gft[4] = {0},
     *gfi[4] = {0};

int open_qti_files()
{
    char fn[1024] = {0};
    int i = 0;
    for (; i < 4; ++i)
    {
        sprintf(fn, "./output/q%d.txt", i);
        gfq[i] = fopen(fn, "w");
        sprintf(fn, "./output/t%d.txt", i);
        gft[i] = fopen(fn, "w");
        sprintf(fn, "./output/i%d.txt", i);
        gfi[i] = fopen(fn, "w");
    }
    return 0;
}

int open_debug_files()
{
    char fn[1024] = {0};
    int i = 0;
    for (; i < 4; ++i)
    {
        sprintf(fn, "./output/fp%d.txt", i);
        gf[i] = fopen(fn, "w");
    }
    return 0;
}

int close_files()
{
    int i = 0;
    for (; i < 4; ++i)
    {
        if (gf[i] != 0)
            fclose(gf[i]);
        if (gfq[i] != 0)
            fclose(gfq[i]);
        if (gft[i] != 0)
            fclose(gft[i]);
        if (gfi[i] != 0)
            fclose(gfi[i]);
    }
    return 0;
}