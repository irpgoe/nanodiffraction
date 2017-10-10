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
 * mex gate_stxm.c stxm.c -output gate_stxm
 */

#include "mex.h"
#include "stxm.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],        /* output variables */
                  int nrhs, const mxArray *prhs[])  /* input variables */
{
    double qrMin;           /* input scalar */
    double qrMax;           /* input scalar */
    double *inData;         /* 1xN input matrix */
    size_t M;               /* size of matrix */
    size_t N;               /* size of matrix */
    double *inMask;         /* 1xN input matrix */
    double *outData;        /* output scalar */

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayStxmWithMask:nrhs","Four inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayStxmWithMask:nlhs","One output required.");
    }
    
    /* get the value of the scalar input  */
    qrMin = mxGetScalar(prhs[2]);
    qrMax = mxGetScalar(prhs[3]);

    /* create a pointer to the real data and mask */
    inData = mxGetPr(prhs[0]);
    inMask = mxGetPr(prhs[1]);

    /* get dimensions of the input data matrix */
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar(0);

    /* get a pointer to the real data in the output matrix */
    outData = mxGetPr(plhs[0]);

    /* call the computational routine */
    if (qrMin != 0 || qrMax != 0){
        stxmWithMask(inData,inMask,outData,(int)M,(int)N);
    }
    else
    {
        stxmWithoutMask(inData,outData,(int)M,(int)N);
    }
}
