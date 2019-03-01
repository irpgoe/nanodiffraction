.. _functions:

##################
External functions
##################

Some generic functionality is contained within seperate functions located inside the ``/functions`` directory. Most functions are used by the ``nanodiffraction.analyze_scan`` routine internally. However, they could also be used stand-alone. Some examples for a stand alone usage of these functions is given further below. Several functions should be mentioned as they are not used by any module but can be nonetheless helpful. ``save_progress`` is a function that is frequently used to store the entire workspace into a single variable and add a comment to a log file for further reference. ``split_struct`` is another helpful routine that searches through structures and nested structures to find a field name that corresponds to a query field name given to the function and returns this value. This function can therefore split values stored in nested structures into a list of values that can be returned as a variable output argument list. 

Examples
========

If we assume that we have only a single detector frame to analyze, it seems as a rather large overhead to invoke ``nanodiffraction.analyze_scan`` for this purpose. It would be clearly much faster and more intuitive to use the pca or azimuthal averaging-routine directly without any larger detour. One example could e.g. look like this:

