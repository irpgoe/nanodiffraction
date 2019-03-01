.. _files:

############
Files module
############


Initialization
==============

The ``files`` module is initialized as in the following example: 

    >>> f = files( 'beamline','id13',...
           'prepath','/path/to/data',...
           'newfile','herz2_roi2',...
           'detector','eiger',...
           'scan',201);
           
By defining ``eiger`` as the type of detector used, the ``eiger`` module is loaded and initialized, see also :ref:`sec_detector_modules`.

Examples
========

Assuming you would like to read the first frame that was acquired for a given file prefix (newfile). To do this, the above section of code would work. Typically, not more than that should be necessary to read data. Once, the files module has been initialized, data can be read using 
    
    >>> data = f.read(1);
    
While this is the most central function that the files class offers, another frequent function that has to be used is ``files.set_<detector>_scan``, e.g. ``files.set_eiger_scan``. These functions define which data set starting with a given prefix (newfile) should be analyzed.

.. _sec_detector_modules:

Detector modules
================

Assuming you would like to read the first frame that was acquired for a given file prefix (newfile). To do this, the above section of code would work. If you would like to refrain from using the ``files`` class, you can more directly use the eiger module instead. This should give you more flexibility if a different task not directly related to scanning experiments should be accomplished instead. For this purpose, you can initialize the eiger module as shown in the following example

    >>> framesPerFile = 126;
    >>> e = eiger(struct('path','/path/to/data') , framesPerFile);
    
For a Pilatus detector, initialization would look something like:

    >>> p = pilatus(struct('filename','my_file_prefix','path','/path/to/data'));
    
Both classes in this example have a method ``read`` that can be used to access the data. The usage is the same as using ``files.read``, e.g.

    >>> data = e.read(1);
    >>> % or
    >>> data = p.read(1);


Parameters
==========

.. list-table:: List of parameters used by the files module
   :widths: 25 75
   :header-rows: 1

   * - Parameter
     - Comment
   * - prepath
     - Path to folder where the data is stored. The prepath depends on the beamline you choose. Check help ``path_by_beamline`` for more information on the default settings. Also check ``show_configuration`` for more information.
   * - newfile
     - File name prefix.
   * - beamline
     - The beamline argument is necessary, because it is referred to when a certain standard file naming convention should be used. See the help on path_by_beamline for more information.
   * - detector
     - Type of detector used in the experiment. Each detector is contained in a separate module, such that each module could also be used independently.
   * - scan
     - Number of the scan. Be aware, that SPEC or other logging software might use a different numbering scheme than the detector. The scan number always refers to the number that is appended to the filename.
   * - fn
     - File number. Note, that some detectors start the numbering scheme with 1, while others start with 0. 
   * - fpf
     - Frames per file. If a container format is used that collects a given number (the fpf) of scattering patterns, this number of frames per container is termed "frames per file". E.g. the Eiger detector uses the hdf5 container format. fpf can be an arbitrary number, but usually it is the number of scan points along the fast axis.
   * - firstFile
     - Defines the first file number that was part of a scan. Although, when the Eiger detector is used this number is usually 1, this is not necessarily the case. The Pilatus detector for example writes single frames instead of blocks and uses a consecutive numbering scheme. Therefore, the first file in the scan has to be indicated. Note that this value is usually set using ``nanodiffraction.set_scan_info``.
   * - slash
     - Used only internally.
   * - detectors
     - Used only internally.
   * - debug
     - 0 by default. If set to 1 can print more information if e.g. a file cannot be read.
   * - implemented
     - Used only internally.


Methods
=======

.. list-table:: List of methods implemented in the files module
   :widths: 25 75
   :header-rows: 1

   * - Method
     - Comment
   * - read
     - Reads a single file. Requires the file number. The file name and location will be determined from the given settings.
   * - set_eiger_scan
     - Routine that updates the eiger class with new scan information.
   * - set_spec_scan
     - Routine that updates the spec class with new scan information.
   * - set_xia_scan
     - Routine that updates the xia class with new scan information.
   * - set_pilatus_scan
     - Routine that updates the pilatus class with new scan information.
   * - which_slash
     - Used only internally.
   * - show_configuration
     - Function that can be used to display the current path settings and all default settings.
   * - remove_leading_trailing_slashes
     - Used only internally.
   * - path_by_beamline
     - Calculates the data path based on the beamline at which the experiment took place. Each beamline usually relies on their own folder structure. This function needs to be extended to accomodate for any uncommon directory structure.
