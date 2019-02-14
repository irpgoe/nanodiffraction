.. _getting_started:

###############
Getting Started
###############


If you want to analyze scattering data, let's get started! 

Installation
============

The toolbox is available from GitHub, see http://irpgoe.github.io/nanodiffraction for more information. Simply clone the toolbox into a directory of your choice using

    >>> git clone https://github.com/irpgoe/nanodiffraction.git

and add the directory to your Matlab search path.


Reading data
============

A few words about the intended structure of the nanodiffraction toolbox. It is logically separated into three main modules, each dedicated to its own task. First, there is the ``files`` module. It can be initialized with the following line of code:

    >>> f = files();
    
Usually, the ``files`` module is used to read a detector image. For this purpose, we use the ``read`` method of the files module, as in the following example:

    >>> f.read(1); % reads the first frame
   
Now, a few questions arise: (1) What detector was used and (2) what does this number (1) stand for? If the files module is initialized as above, default values are used for the path to the location where the data is stored, the type of detector that was used, etc. Let us therefore look at a more realistic example:
 
    >>> f = files( 'beamline','id13',...
           'prepath','/home/Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2018/extern/ESRF_ID13_ls2728_Jan2018/DATA/AUTO-TRANSFER/eiger1',...
           'newfile','setup3',...
           'detector','eiger',...
           'scan',11);   
           
This set of values is the essential minimum that the files module need to access the data. A complete description is given in the respective chapter on the files module.
Now, you should be able to read a single detector frame using

    >>> frame = f.read(1); % reads the first frame
    
The index, here 1, is the frame number, i.e. the n-th (here: first) frame acquired for a given ``newfile``.

Processing data
===============

To actually analyze the data, we need to have methods that can be applied to process the data obtained from ``files``. We use the ``nanodiffraction`` module for this purpose. It is basically a collection of all geometrical and experimental parameters and related methods that can be used for batch processing of the data. This module can be either initialized with default values or directly with the appropriate experimental parameters, as in the following lines of code:

    >>> e = nanodiffraction('energy',13E3,...
           'detDistance',0.971434,...
           'pixelsize',75E-6,...
           'Ny',2070,...
           'Nz',2167,...
           'pby',1468.934,...
           'pbz',1315.318);  
           
The parameters should in this case be rather self-explanatory. For a more in-depth explanation of each parameter, see the respective module description. Let us proceed to analyze a simple STXM (scanning transmission X-ray microscopy) scan. For this purpose, we want to go through all frames within a single scan and calculate the total number of scattered photons. First, we need to initialize the scan parameters with the following line:

    >>> e.set_scan_info('SNy',100,'SNz',100);
    
This defines the horizontal (SNy) and vertical (SNz) scan dimensions. We then simply need to call the most used method of this module, ``analyze_scan``:

    >>> e.analyze_scan('method','stxm');
    
But, hugh, it doesn't work? That's because the ``nanodiffraction`` module does not yet know, where to get the data from. The toolbox therefore contains a function that can link the processing module to the file reading module:

    >>> link(f,e)
    
A full (but minimal) example might therefore look like this: 

    >>> f = files(...);
    >>> e = nanodiffraction(...);  
    >>> link(f,e)
    >>>    
    >>> e.set_scan_info('SNy',100,'SNz',100);
    >>> result = e.analyze_scan('method','stxm');    
    
    
Visualizing data
================

For simple visualization purposes, you might use built-in Matlab functions (imagesc, imshow, etc...). You could however also make use of the ``display`` module that contains pre-defined routines for visualization of diffraction data and resulting 2D maps of some structure parameter.

In general, it is as simple as typing:

    >>> d = display();
    >>> d.stxm(darkfield);
    
    
Summary
=======

The above is short description of how a simple use case of the nanodiffraction toolbox could look like. More information is given in the Guide section. In addition, each function contains a help block. The usage of each function can therefore be obtained in Matlab by typing

    >>> help pca
    
for the case of the function ``pca``. You can also obtain help on a module

    >>> help nanodiffraction
    
or a method of a module

    >>> help files.show_configuration
    
In all other cases, feel free to contact the author of the toolbox directly.


    