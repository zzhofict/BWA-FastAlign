/*********************************************************************************************
    Description: data and files for debugging
    Copyright : All right reserved by .

    Author : 
    Date : 2024/04/09
***********************************************************************************************/
#include <stdio.h>

////////////////// for debug and test //////////////////////////

#define DEBUG_FILE_OUTPUT // gfp1-4，debug
// #define COUNT_SEED_LENGTH // seed1，
// #define GET_FULL_MATCH_READ // reads
// #define COUNT_CALC_NUM // BSW
// #define GET_DIFFERENT_EXTENSION_LENGTH // extensionquery，target，
// #define GET_KSW_ALIGN_SEQ
// #define DEBUG_SW_EXTEND // bswdebug

////////////////////////////////////////////////////////////////

extern FILE *gf[4], *gfq[4], *gft[4], *gfi[4];

int open_qti_files();

int open_debug_files();

int close_files();