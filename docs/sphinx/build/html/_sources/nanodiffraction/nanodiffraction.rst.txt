.. _nanodiffraction:

######################
Nanodiffraction module
######################

Initialization
==============

The ``nanodiffraction`` module is initialized as in the following example: 

    >>> e = nanodiffraction('energy',14.85E3,...
                    'detDistance',0.971434,...
                    'pixelsize',75E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1468.934,...
                    'pbz',1315.318);


Examples
========


Parameters
==========

.. list-table:: List of parameters used by the nanodiffraction module
   :widths: 25 75
   :header-rows: 1

   * - Parameter
     - Comment
   * - data
     - Contains a handle to a ``files`` class. This handle is set using the function ``link`` or could be set manually.
   * - helper
     - Structure containing short helper functions that have proven useful in analyzing data and as shorthand notations.
   * - scan
     - Structure containing scan parameters. The values in this structure are defined using ``nanodiffraction.set_scan_info``.
   * - energy
     - Energy of the X-radiation in eV.
   * - wavelength
     - Wavelength  of the X-radiation in m.
   * - wavenumber
     - Wavenumber of the X-radiation in m⁻¹.
   * - yaxis
     - Detector y-axis in pixel units.
   * - zaxis
     - Detector z-axis in pixel units.
   * - R
     - 2D array storing the pixel-wise distance to the primary beam in pixel units.
   * - theta
     - 2D array storing the pixel-wise scattering angle with respect to the primary beam in radians.
   * - Ny
     - Number of pixels along y coordinate.
   * - Nz
     - Number of pixels along z coordinate.
   * - pby
     - Primary beam position along y coordinate in pixels.
   * - pbz
     - Primary beam position along z coordinate in pixels.
   * - pixelsize
     - Size of a single (square) pixel in m.
   * - detDistance
     - Focus to detector distance.
   * - qyaxis
     - y-axis of the detector in units of m⁻¹.
   * - qzaxis
     - z-axis of the detector in units of m⁻¹.
   * - qr
     - 2D array storing the pixel-wise distance to the primary beam in units of m⁻¹.
   * - phi
     - 2D array storing the pixel-wise azimuthal angle with respect to the -y coordinate in degrees [0 to 360 degrees].
   * - phi2
     - 2D array storing the pixel-wise azimuthal angle with respect to the -y coordinate [0 to 180 degrees].
   * - Ny_orig
     - Number of pixels along y coordinate of the raw data (without binning or cropping of the detector).
   * - Nz_orig
     - Number of pixels along z coordinate of the raw data (without binning or cropping of the detector).
   * - pby_orig
     - Primary beam position along y coordinate in pixels (without binning or cropping of the detector).
   * - pbz_orig
     - Primary beam position along z coordinate in pixels (without binning or cropping of the detector).
   * - pixelsize_orig
     - Size of a single (square) pixel in m (without binning or cropping of the detector).
   * - mask
     - 2D logical array defining bad pixels in the detector images. A (1) corresponds to an invalid detector pixel, while a (0) corresponds to a valid detector pixel. The mask is set using ``nanodiffraction.set_mask``.
   * - sels
     - 3D logical array defining valid pixels that should be taken into account for the anlysis. Invalid pixels (0) are not analyzed. A selection is initially defined using ``nanodiffraction.set_selection``. More selections can be added using ``ǹanodiffraction.add_selection``.
   * - corr
     - 2D array that can be multiplied with the original data frame to correct for the absorption by semi-transparent objects such as beamstop holders. The correction matrix is set using ``nanodiffraction.set_correction``.
   * - mask_orig
     - Same as mask, however, with the size of the original data frame (without binning or cropping of the detector).
   * - sels_orig
     - Same as sels, however, with the size of the original data frame (without binning or cropping of the detector).
   * - corr_orig
     - Same as corr, however, with the size of the original data frame (without binning or cropping of the detector).
   * - empty
     - Detector image that should only contain background signal that can be subtracted from the raw data upon batch analysis. The empty image is set using ``nanodiffraction.set_empty``.
   * - detectorRoi
     - Can be either 'on' or 'off'. If set to 'on', every data frame will be automatically cropped to a region defined by 'detRoiY' and 'detRoiY'. 
   * - binning
     - Can be either 'on' or 'off'. If set to 'on', every data frame will be automatically binned using the binning ratio defined by 'biny' and 'binz'.
   * - biny
     - Integer binning ratio along y-axis.
   * - binz
     - Integer binning ratio along z-axis.
   * - detRoiY
     - Two element vector containing the cropping limits along the y-axis. The largest possible value is [1 Ny].
   * - detRoiZ
     - Two element vector containing the cropping limits along the z-axis. The largest possible value is [1 Nz].
   * - fluoCounters
     - Cell array of counters that define regions of interests for certain elements in a 1D fluorescence detector. In units of channel numbers. E.g. ``fluoCounters = {{'total',100,2600},{'feka',450,500}};``
     
Methods
=======

.. list-table:: List of methods implemented in the nanodiffraction module
   :widths: 25 75
   :header-rows: 1

   * - Method
     - Comment
   * - add_selection
     - Can be used once ``set_selection`` was used to initially define a selection. Adds a selection to the set of selections. In batch processing, all selections are sequentially applied on the same data frame such that different regions on the detector can be analyzed at once, without reloading of data.
   * - analyze_scan
     - Performs batch processing of a scan with scan dimensions set by ``nanodiffraction.set_scan_info``. Loops through all data frames, performs the required analysis routines for each selection separately and returns a structured data output. 
   * - around_y
     - Useful shorthand routine when e.g. setting a detector roi. Returns a two-element vector with entries [pby-X pby+X], where X is an integer value in pixel units. 
   * - around_z
     - Useful shorthand routine when e.g. setting a detector roi. Returns a two-element vector with entries [pbz-X pbz+X], where X is an integer value in pixel units.
   * - azimuthal_mask
     - Defines an azimuthal mask with start and end values 'phi1' and 'phi2' in degrees.
   * - azimuthal_spread
     - Interpolates a one-dimensional curve onto the detector by azimuthal spreading.
   * - calculate_composite
     - Calculates a composite image. Usually visualized in conjunction with ``display.composite``
   * - calculate_empty
     - Wrapper for ``nanodiffraction.analyze_scan`` that can be used to calculate an average image. This image can be used as an empty image within ``nanodiffraction.analyze_scan`` using ``nanodiffraction.set_empty``.
   * - clicktool
     - One a STXM map is visualized in the active figure, clicktool can be executed which activates a mouse courses such that one can click into the 2D map. By clicking onto a scan point in the map, a diffraction image will be displayed that was collected at this location. 
   * - copy
     - The nanodiffraction class can be copied. As an example, one could create a template class ``templ = nanodiffraction(...)`` and from this template, a SAXS and a WAXS class can be obtained ``saxs = copy(templ)`` or ``waxs = copy(templ)``.
   * - get_index
     - 
   * - get_location_in_frame
     - 
   * - get_parameters
     -  
   * - heal
     - 
   * - ind2sub
     - 
   * - make_mask_symmetric
     - 
   * - moment_from_selection
     - 
   * - n_of_q
     - Converts q coordinates into pixel coordinates.
   * - process
     -  
   * - process_mask
     - 
   * - q_of_n
     - Converts pixel coordinates into q coordinates.
   * - qr_tilt
     - 
   * - radial_mask
     - 
   * - read
     - 
   * - readm
     -            
   * - readp
     -       
   * - refresh
     -                       
   * - remove_semitransparent_object
     -            
   * - set_corr
     -       
   * - set_empty
     -       
   * - set_mask
     -            
   * - set_roi_and_binning
     -       
   * - set_scan_info
     -       
   * - set_selection
     -            
   * - simulate_actomyosin
     -       
   * - simulate_agbh
     -       
   * - simulate_airscattering
     -            
   * - simulate_bs
     -       
   * - simulate_streak
     -       
   * - stitch
     -            
   * - strip
     -       
   * - sub2ind
     -       
   * - verify_array
     -            
                                 
     