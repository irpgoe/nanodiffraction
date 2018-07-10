function [data,low,high] = autoc(data)
% AUTOC  estimates the contrast of an image to lie within the 5 and 95
% percentile of an image. Data will be clipped to the according value
%
%   [DATA,LOW,HIGH] = AUTOC(DATA)
%
%   The following arguments are required:
%       DAT:: ()
%           Diffraction pattern.
%
%       LOW:: ()
%           lower percentile.
%
%       HIGH:: ()
%           upper percentile.
%
% Copyright 2017 Institute for X-ray Physics (University of Göttingen)

% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
% OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    low  = prctile(reshape(data,1,[]),5);
    high = prctile(reshape(data,1,[]),95);

    % clipping
    data(data < low) = low;
    data(data > high) = high;
end