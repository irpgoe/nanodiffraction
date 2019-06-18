function image = binning_2d(image,binx,biny)                
% BINNING_2D  bins and image according to the binning
% factors binx and biny.
%
%   ``BINNED_IMAGE = BINNING_2D(IMAGE, BINX, BINY)``
%
% Parameters
% ----------
%   image: Two-dimensional numeric array
%       Input image
%
%   binx: Integer
%       Integer binning factor along the horizontal dimension
%
%   biny: Integer
%       Integer binning factor along the vertical dimension
%
% Returns
% -------
%   binned_image: Two-dimensional numeric array
%       Image binned by the binning factors binx and biny
%
% Notes
% -----
% Example 1:   
%
% .. code-block:: matlab
%
%   testimage = magic(256);
%   binned_im = binning_2d(testimage,2,2);
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

    [m,n] = size(image);
    nnew = n - mod(n,binx);
    mnew = m - mod(m,biny);

    % now it is divisible by binx, biny
    image = image(1:mnew,1:nnew);

    % binning along z (correct)
    image=sum( reshape(image,biny,[]),1); 
    image=reshape(image,mnew/biny,[]); 

    % binning along y
    image=sum(reshape(image',binx,[]),1); 
    image=reshape(image,nnew/binx,[])';
    
end