/*==========================================================
 * pca.c - Principal component analysis of diffraction data
 *
 * Calculates the two main components of the data
 * by Principal Component Analysis, as detailed in
 * Bernhardt et.al (2016), Biophys. J.
 *
 * The calling syntax is:
 *
 *		pca(data, mask, qx, qy, w, phi, m, n)
 *
 * Copyright 2017 Institute for X-ray Physics (University of Göttingen)
 *
 * Permission is hereby granted, free of charge, to any person obtaining 
 * a copy of this software and associated documentation files (the "Software"), 
 * to deal in the Software without restriction, including without limitation 
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 * and/or sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *========================================================*/

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define PI 3.14159265

/* The computational routine */
void pca(double *x, double *mask, double *qx, double *qy, double *w, double *phi, int m, int n)
{    
    int col;
    int row;
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
