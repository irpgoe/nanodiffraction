.. _mex:

#############
Mex-functions
#############

Mex-functions are pre-compiled C-functions that can significantly improve computational speed. Mex-functions are used primarily in the context of batch processing where a small increase in a single computation, e.g. the computation of the orientation of the scattering signal, can significantly impact overall computation time. At this time, there are only a small number of mex functions available, see the following list:

    * STXM
    * PCA
    
Mex functions were compiled on a linux machine. Below is an example where we have compiled the function pca.c implementing a principal component analysis for the special case of scattering patterns. The C-function was compiled using the ``gcc```compiler and the following command:

    >>> gcc -Wall -fPIC -c pca.c -lgsl -lgslcblas -o pca.o
    
The compiled functions were the packaged into a static library using

    >>> ar -cvq libpca.a pca.o
    
You can check the contents of the built library by using ar: 

    >>> ar -t libpca.a
    
At last, we need to compile the mex functions. This can be done using the command ``mex``:
   
   >>> mex ./toolbox/nanodiffraction/mex/mex_pca_test.c -Ltoolbox/nanodiffraction/mex -lpca -lgsl -lgslcblas
   
One common error when executing the mex command is, that the arguments were passed in the wrong order. Since e.g. pca.c relies on the gsl library, when -lpca is called to resolve the function ``pca``, it also tries to resolve libgsl objects. If they were included before -lpca, these objects stay unresolved. 

Pre-compiled mex-functions are available from within the toolbox and should be located in the ``/mex`` directory.