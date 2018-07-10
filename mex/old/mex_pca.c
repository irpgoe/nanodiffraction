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
 * mex -lgsl -lgslcblas ./toolbox/nanodiffraction/mex/mex_pca.c 
 */

#include "mex.h"
/* #include <stdio.h> */
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define PI 3.14159265

/* The computational routine */
void pca(double *x, double *mask, double *qx, double *qy, double *w, double *phi, mwSize m, mwSize n)
{
    /*
    % first order moments
    tr = sum(dat(:));
    
    % calculate covariance
    cov_qx_qx= sum(sum(dat.*qy.^2))./tr;
    cov_mix  = sum(sum(dat.*qy.*qz))./tr; % Mischterm
    cov_qy_qy= sum(sum(dat.*qz.^2))./tr;

    % covariance matrix
    A=[cov_qx_qx cov_mix; cov_mix cov_qy_qy];

    % solve eigenvalue problem (V: Eigenvectoren, D:
    % Diagonalmatrix)
    [V,D] = eig(A);
    */
    
    mwSize col;
    mwSize row;

    double tr = 0;
    
    /* total sum of x */
    for (row=0; row<m; row++) {
        for (col=0; col<n; col++) {
            x[row*m + col] = x[row*m + col] * mask[row*m + col];
            tr += x[row*m + col];
        }
    }
    
    /* covariance matrix */
    double cov_qx = 0;
    double cov_mixed = 0;
    double cov_qy = 0;
    
    /* total sum of x */
    for (row=0; row<m; row++) {
        for (col=0; col<n; col++) {
            cov_qx += x[row*m + col] * pow(qx[row*m + col],2);
            cov_qy += x[row*m + col] * pow(qy[row*m + col],2);
            cov_mixed += x[row*m + col] * qx[row*m + col] * qy[row*m + col];
        }
    }
    
    double data[] = {cov_qx/tr   , cov_mixed/tr,
                     cov_mixed/tr, cov_qy/tr   };
    
    gsl_matrix_view mv 
        = gsl_matrix_view_array (data, 2, 2);

    gsl_vector *eval = gsl_vector_alloc (2);
    gsl_matrix *evec = gsl_matrix_alloc (2, 2);     
    
    gsl_eigen_symmv_workspace * ws
        = gsl_eigen_symmv_alloc (2);
  
    gsl_eigen_symmv (&mv.matrix, eval, evec, ws);

    gsl_eigen_symmv_free (ws);

    gsl_eigen_symmv_sort (eval, evec, 
                        GSL_EIGEN_SORT_ABS_DESC);
  
    /* get eigenvalues and eigenvectors (assuming they are sorted) */
    double eval_i = gsl_vector_get (eval, 0);
    double eval_j = gsl_vector_get (eval, 1);
       
    gsl_vector_view evec_i 
       = gsl_matrix_column (evec, 0);
    
    double xcomp = gsl_vector_get (&evec_i.vector, 0);
    double ycomp = gsl_vector_get (&evec_i.vector, 1);
      
    /* anisotropy */
    *w = (eval_i - eval_j) / (eval_i + eval_j);
    
    /* orientation */
    *phi = atan2(ycomp,xcomp) * 180.0 / PI; 
    if (*phi < 0){
        *phi += 180;
    }
    
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    
    /* debugging information */
    /*
    printf ("%g\t%g\t%g\n", cov_qx/tr, cov_mixed/tr, cov_qy/tr );

    printf ("eigenvalue 1 = %g\n", eval_i);
    printf ("eigenvalue 2 = %g\n", eval_j);
        
    printf ("eigenvector = \n");
    gsl_vector_fprintf (stdout, 
                        &evec_i.vector, "%g");
    
    printf ("x = %g\n", xcomp);
    printf ("y = %g\n", ycomp);    
    printf ("w = %g\n", *w);
    printf ("phi = %g\n", *phi);
    */
}

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
    pca(inData,inMask,qx,qy,w,phi,(mwSize)M,(mwSize)N);

}
