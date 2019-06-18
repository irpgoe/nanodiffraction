function [x,sortIndex,n,bin] = b1d_init(grid,mask,sel,bins)
% B1D_INIT  calculates the mapping between detector pixels and q-bins.
%   
%   ``result = b1d_init(grid,mask,selection,bins)``
%
% Parameters
% ----------
%   grid: Two-dimensional numeric array
%       Detector pixel coordinates
%
%   mask: logical array, default = []
%       Logical array of the size of the detector grid that defines which pixels are considered bad (1 = bad, 0 = valid)
%
%   selection: logical array, default = []
%       Logical array of the size of the detector grid that defines which pixels should be analyzed (1 = valid, 0 = invalid)
%
%   bins: Numeric value, default = 360, optional
%       Number of bins used for histogramming
%
% Returns
% -------
%   x: One-dimensional numeric array
%       Contains the x-coordinate based on the grid used
%
%   sortIndex: One-dimensional numeric array
%       Contains the indices of the q-bins for each detector pixel
%
%   n: One-dimensional numeric array
%       Contains the number of pixels that belong to each q-bin
%
%   bin: One-dimensional numeric array
%       bin index
%
% Notes
% -----
% Example 1:
%
% .. code-block:: matlab
%
%   e = nanodiffraction();
%   testdata = e.qr;
%   [x,sortIndex,n,bin] = b1d_init(grid,mask,sel,bins)
%   testavg = b1d(testdata,[],[],sortIndex,n,bin);
%
% See also B1D, B1D_FAST, ANALYZE_SCAN

    if isempty(bins)
        bins = 360;
    end

    % sort grid
    x = grid(~(mask | ~sel)); 
    [x, sortIndex] = sort(x);
    
    % binning of grid
    binEdge = linspace(min(x),max(x),bins);
    [n,~,bin] = histcounts(x,binEdge);
    
    % remove trailing zeros
    x = binEdge(1:sum(n~=0));
    
end