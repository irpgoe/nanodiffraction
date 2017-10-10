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
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#define PI 3.14159265

/* The computational routine */
void symmetry(double *x, double *bgr, double sigma, double noise_level, int m, int n)
{    
    int delta_theta = 20; /* hard coded theta window size */
    
    int col;
    int row;
        
    /* calculate snr */
    for (row=0; row<m; row++) {
        for (col=0; col<n; col++) {
            if (bgr[row*m + col] == 0){
                x[row*m + col] = 0;
            } else {
                x[row*m + col] = (x[row*m + col] - bgr[row*m + col]) / sqrt(bgr[row*m + col]);
        }
    }
    
    /* next we implement 
% Gaussian filter (remove noise)
h = fspecial('gaussian', [10 10], sigma);
h = h./sum(h(:));
snr = conv2(snr,h,'same');

% binary image (snr threshold)
binary_image = (snr >= noise_level);
        */
        
    /* prepare fft */
    work = gsl_fft_real_workspace_alloc (n);
    real = gsl_fft_real_wavetable_alloc (n);
        
    /* fourier transform the rows */
    for (int row = 0; row < m; row++){
        double* data = x[row*m];
        gsl_fft_real_transform (data, 1, n, real, work);
        gsl_fft_real_wavetable_free (real);
        /* apply filter HERE */
        hc = gsl_fft_halfcomplex_wavetable_alloc (n);
        gsl_fft_halfcomplex_inverse (data, 1, n, hc, work);
        gsl_fft_halfcomplex_wavetable_free (hc);
        gsl_fft_real_workspace_free (work);
    }
    /* fourier transform the rows */
        
	/* apply gauss filter */
        
	/* fourier-backtransform the rows */
    /* fourier-backtransform the cols */
        
	/* apply threshold */
    for (row=0; row<m; row++) {
        for (col=0; col<n; col++) {
            if (x[row*m + col] > noise_level){
                x[row*m + col] = 1;
            } else {
                x[row*m + col] = 0;
        }
    }
        
        
}
