.. _display:

##############
Display module
##############

Initialization
==============

The ``display`` module is initialized as in the following example: 

    >>> d = display();
           
           
Examples
========

Assuming you have a scattering map, e.g. darkfield contrast, and want to display the map with a proper y- and z-scale and a default STXM colormap. In this case, the following should be sufficient: 
    
    >>> d.stxm(darkfield_map, struct('sampl',10,'scale',1,'units','um');
    
This places y- and z-ticks at every 10th. scan point. The step size was 1 µm, given by the ``scale`` and ``units`` parameter. Now, if you prefer a scale bar an no axis labels (like me), you make adapt the example slightly:

    >>> d.stxm(darkfield_map, struct('sampl',10,'scale',1,'units','um');
    >>> d.remove_axis();
    >>> d.add_scalebar(10); 
    >>> d.round_limits(3);
    
This adds a scale bar with length of 10 steps. In this case, this corresponds to 10 µm. The function ``round_limits`` is used to get nice, round limits of the color axis.  


Parameters
==========

.. list-table:: List of parameters used by the display module
   :widths: 25 75
   :header-rows: 1

   * - Parameter
     - Comment
   * - 
     - not yet complete.
   

Methods
=======

.. list-table:: List of methods implemented in the files module
   :widths: 25 75
   :header-rows: 1

   * - Method
     - Comment
   * - 
     - not yet complete.