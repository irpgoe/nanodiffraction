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

/* Compiled with
 * mex -lgsl -lgslcblas gate_pca.c pca.c -output gate_pca
 */

#include "mex.h"
/* #include <stdio.h> */
#include "pca.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],        /* output variables */
                  int nrhs, const mxArray *prhs[])  /* input variables */
{
    /* double placeholder; */   /* input scalar */ 
    double *inData;         /* 1xN input matrix */
    double *inMask;         /* 1xN input matrix */
    double *qx;             /* 1xN input matrix */
    double *qy;             /* 1xN input matrix */
    size_t M;               /* size of matrix */
    size_t N;               /* size of matrix */
    double *w;              /* output scalar or matrix */
    double *phi;            /* output scalar or matrix */

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayStxmWithMask:nrhs","Four inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayStxmWithMask:nlhs","Two outputs required.");
    }
    
    /* get the value of the scalar input  */
    /* qrMin = mxGetScalar(prhs[1]); */

    /* create a pointer to the data */
    inData  = mxGetPr(prhs[0]);
    inMask  = mxGetPr(prhs[1]);
    qx      = mxGetPr(prhs[2]);
    qy      = mxGetPr(prhs[3]);

    /* get dimensions of the input data matrix */
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar(0);
    plhs[1] = mxCreateDoubleScalar(0);

    /* get a pointer to the real data in the output matrix */
    w   = mxGetPr(plhs[0]);
    phi = mxGetPr(plhs[1]);

    /* call the computational routine */
    pca(inData,inMask,qx,qy,w,phi,(int)M,(int)N);

}
