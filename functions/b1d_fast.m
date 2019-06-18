function [results] = b1d_fast(dat,mask,sel,sortIndex,n,bin)
% B1D_FAST  performs an angular average based on histogramming the intensity of
% each pixel in a certain radial interval. b1d is short for one-dimensional binning and
% is considered to be a rather simple but fast method to get accurate
% angular averages. B1D_FAST is used in batch processing of files and based
% on the method B1D, however, it avoids the repeated calculation of how
% many detector pixels belong into one q-bin. Instead, this needs to be
% calculated beforehand and passed to B1D_FAST.
%   
%   ``result = b1d_fast(data,mask,selection,sortIndex,n,bin)``
%
% Parameters
% ----------
%   data: Two-dimensional numeric array
%       Diffraction pattern
%
%   mask: logical array, default = []
%       Logical array of the size of the input data that defines which pixels are considered bad (1 = bad, 0 = valid)
%
%   selection: logical array, default = []
%       Logical array of the size of the input data that defines which pixels should be analyzed (1 = valid, 0 = invalid)
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
% Returns
% -------
%   result: structure
%       Structure with the following fields:
%
%       - dat_1d: One-dimensional azimuthally integrated data
%       - y: A short-hand for dat_1d
%       - error: Propagated measurement error, calculated from the square
%       root of the average intensity I, divided by the square root of the
%       number of data points n averaged in each bin: sqrt(I)/sqrt(n)
%       [point-wise division]
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
% See also B1D, B1D_INIT, ANALYZE_SCAN

% Copyright 2017 Institute for X-ray Physics (University of Göttingen)
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
% associated documentation files (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge, publish, distribute,
% sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or
% substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
% NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
% OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
            
    % Settings, if none are passed than default values are taken
    if isempty(sel)
        sel = ones(size(dat));
    end
    if isempty(mask)
        mask = zeros(size(dat));
    end
    
    % put intensity values into bins (each qr value has a bin associated) and take mean
    aavg = dat(~(mask | ~sel));
    aavg = aavg(sortIndex);
    aavg = accumarray(bin,aavg,[],@mean,NaN);
    
    % remove NaNs from intensity and bin index (NaNs appear where n is 0)
    n(isnan(aavg)) = [];
    aavg(isnan(aavg)) = [];
    
    % save results
    results.dat_1d = aavg;
    results.y = aavg;
    results.error = sqrt(abs(aavg))./sqrt(n');
end