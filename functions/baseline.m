function data = baseline(data, winSize, iterations)
% BASELINE  Iterative baseline identification using the strip algorithm as used in pyMCA, see Kieffer et al., 2014.
%   
%   ``data = baseline(data, winSize, iterations)``
%
% Parameters
% ----------
%   data: One-dimensional numeric array
%       Input data
%
%   winSize: integer
%       Integer value greater than 0, that defines the size of the filter window used for the strip algorithm.
%
%   iterations: integer
%       Integer value greater than 0, that defines the number of iterations. 
%
%
% Returns
% -------
%   data: One-dimensional numeric array
%       Baseline that can be used for background subtraction
%
% Notes
% -----
% Example 1:
%
% .. code-block:: matlab
%
%       [x,y]   = split_struct(b1d(data,[],[],e.qr,512),{'x','y'});
%       bl      = baseline(y, 4, 1000);
%
% Example 2:
%
% .. code-block:: matlab
%
%       e       = nanodiffraction();
%       [x,y]   = split_struct(b1d(data,[],[],e.qr,512),{'x','y'});
%       bl      = baseline(y, 4, 1000);
%       bgr_2d  = e.azimuthal_spread(...);
%       bgr_sub = data - bgr_2d;
%

% Copyright 2017 Institute for X-ray Physics (University of Göttingen)
%
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
% OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
    
    m = data;
    
    for ii = (1+winSize):(numel(data)-winSize)
        m(ii) = 0.5*(data(ii-winSize) + data(ii+winSize));
    end
    sel = m < data;
    data(sel) = m(sel);

    if iterations == 0
        return;
    end
    data = baseline(data, winSize, iterations-1);
end