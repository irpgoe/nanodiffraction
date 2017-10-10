/*==========================================================
 * arrayStxmWithMask.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayStxmWithMask(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "stxm.h"

/* The computational routine */
void stxmWithMask(double *x, double *y, double *out, int m, int n)
{
    int col;
    int row;
    
    /* multiply each element y by x */
    for (row=0; row<m; row++) {
        for (col=0; col<n; col++) {
            *out += x[row*m + col] * y[row*m + col];
        }
    }
}

/* The other computational routine */
void stxmWithoutMask(double *x, double *out, int m, int n)
{
    int col;
    int row;
    
    /* multiply each element y by x */
    for (row=0; row<m; row++) {
        for (col=0; col<n; col++) {
            *out += x[row*m + col];
        }
    }
}
