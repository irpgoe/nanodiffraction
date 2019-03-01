.. _guide:

#####
Guide
#####

In this guide we will analyze data available online [Nicolas2017b]_ and described in [Nicolas2017a]_. Once the data was downloaded, we can start Matlab. At this stage, there is however already one important caveat. The available h5 files are compressed using LZ4 compression. This filter is not part of the standard filter set available from the h5 library, see https://support.hdfgroup.org/services/filters.html. The filter is however available for download. On a linux machine, once the filter is downloaded, it has to be made available by exporting the path to the filter  

    >>> HDF5_PLUGIN_PATH = /path/to/plugin
    
Now, Matlab can be started as usual. Unfortunately, this procedure has only been tested for linux-based operating systems. 

Before proceeding any further, we add the nanodiffraction toolbox to the Matlab search path.

    >>> addpath(genpath('/path/to/nanodiffraction'));

Let us now see, how the ``files`` module predefines certain folder structures for a given set of beamlines

    >>> f.files();
    >>> f.show_configuration(struct('allDefaults',1));
    
Since it is reported that data was collected at beamline ``id13`` using an Eiger detector (Dectris, Switzerland), we have a closer look at the following line:

    >>> ...
    >>> eiger		||	id13		||	/<testprepath>/<testnewfile>_1_master.h5
    >>> ...
    
Then, we initialize the ``files`` module using the correct ``prepath`` argument, as suggested by the above line

    >>> f = files( 'beamline','id13',...
           'prepath','/path/to/data',...
           'newfile','herz2_roi2',...
           'detector','eiger',...
           'scan',201);   
           
and the ``nanodiffraction`` module

    >>> e = nanodiffraction('energy',14.6407E3,...
           'detDistance',1.928650,...
           'pixelsize',75E-6,...
           'Ny',2070,...
           'Nz',2167,...
           'pby',1500.499,...
           'pbz',1372.509);     
                    
and finally the ``display`` module

    >>> d = display();

Now, as already mentioned in the Getting Started section, we have to `tell the nanodiffraction module where the data is stored` using

    >>> link(f,e);
    
as `tell the visualization module what the scanning geometry was` using
    
    >>> link(e,d)
    
Both lines can actually be combined into one:

    >>> link(f,e,e,d)

Because no detector is perfect and contains some sort of artifacts, we have to mask out pixels with unwanted values. In this example, the detector contains modular gaps and hot pixels that we would like to mask out. Fortunately, this is easy to accomplish

    >>> e.set_mask(f.read(1)==(2^32-1)); % read first frame, identify all hot pixels and set this logical array as a mask
    
The function ``set_mask()`` is used to store the detector mask for any later processing. It can always be retrieved using

    >>> current_mask = e.mask;
    
Following data masking we not only want to throw away unwanted pixels, we also would like to highlight which set of pixels should actually be analyzed. In diffraction, this is usually an azimuthal wedge or a radial section. In this example, we use a radial section. Only pixels between 86 and 230 (in radial pixel units) should be analyzed:
    
    >>> e.set_selection(e.radial_mask('r1',86,'r2',230));
    
Here, we actually combined to important functions. The radial mask is defined using ``radial_mask()``. ``r1`` and ``r2`` are in pixel units. If any other unit should be used, it is useful to have a look at the help function to understand the available options, see e.g.
   
    >>> help e.radial_mask
    >>> % help nanodiffraction.radial_mask also works
    
At last, we want to bin data and only analyze data in a rectangular region-of-interest. This is done with the following line:

    >>> e.set_roi_and_binning('binning','on',...
        'biny',4,'binz',4,...
        'detectorRoi','on',...
        'detRoiY',around_y(200),...
        'detRoiZ',around_z(200));
        
Binning and ROI's slightly increase the complexity of the analysis. This is because the nanodiffraction class now actually needs to take care of two geometries, the original geometry and the `virtual` binned and cropped detector. Inspecting the current module

    >>> e
    
You can find e.g. the following lines:
   ...
   pby: 1.4689e+03
   pbz: 1.3153e+03
   ...               
   pby_orig: 1.4689e+03
   pbz_orig: 1.3153e+03
   ...
   
In general, this should not pose any difficulties. Whenever the function ``set_roi_and_binning`` is called, it automatically updates the geometry. Therefore, try the following line

    >>> e.set_roi_and_binning('binning','on',...
        'biny',2,'binz',2,...
        'detectorRoi','on',...
        'detRoiY',around_y(100),...
        'detRoiZ',around_z(100));
    >>> e
    
Did you see how the geometry was updated?

Now, we can start analyzing the data. The last step before actually doing so requires to define the dimensions of the scan:
   
    >>> e.set_scan_info('SNy',101,'SNz',101);
    
The data can now be processed using the ``analyze_scan`` function, the most central piece of the nanodiffraction toolbox:    
    
    %% analyze scan / methods can be combined using the '+' notation
    % as methods: stxm | pca | crystal | symmetry | average | heal | sum can be
    % used
    >>> result = e.analyze_scan('method','stxm+pca');

The output of this function is structured and explained in the corresponding help section (``help e.analyze_scan``). The toolbox offers a helper function that can simplify the structure of the output slightly by splitting it up into separate variables, e.g.

    >>> [df,angle,w] = split_struct(result,{'df','angle','w'});
    
Now we can visualize the data using 
    >>> close all
    >>> f1 = figure(1);
    >>> d.stxm(df,struct('scale',0.02,'sampl',10,'unit','mm','alpha',df>2e3));
    >>> caxis([5000 10000]);
    >>>
    >>> f2 = figure(2);
    >>> d.stxm(angle,struct('scale',0.02,'sampl',10,'unit','mm','alpha',df>2e3));
    >>> hold on;
    >>> d.pca(angle);
    >>>
    >>> f3 = figure(3);
    >>> d.stxm(angle,struct('scale',0.02,'sampl',10,'unit','mm','alpha',df>2e3));
    >>> hold on;
    >>> d.pca(res.pca.angle,struct('quiver','on'));

In fact, the following line from the above block shows how almost all toolbox functions are structured:

    >>> d.stxm(df,struct('scale',0.02,'sampl',10,'unit','mm','alpha',df>2e3));
    
All required arguments are passed directly as arguments to the function. It would therefore also suffice to use

    >>> d.stxm(df)
    
which is in this case the bare minimum for a first output. To tune the output to the actual scanning dimensions and to mask out values that are outside the sample, additional arguments can be passed to the funtion in a single structure, that usually is pre-defined with default values. Note, that not all default values have to be overwritten. See the Matlab help (``help d.stxm``) for more information.

