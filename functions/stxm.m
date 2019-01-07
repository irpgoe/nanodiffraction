function result = stxm(data,mask,sel)
% STXM  Sums a 2d array. In addition, it takes into
% account a minimum and maximum q_r value for integration
%
%   result = stxm(data,mask,selection)
%
% The following arguments are supported:
%       data:: [] (required) 
%           A two-dimensional diffraction pattern. 
%
%       mask:: [] (required)
%           Mask that defines which pixels are bad and should not be taken
%           into account. The mask can be an empty matrix. In this case, it
%           is assumed that all pixels are valid. The dimensions of mask
%           have to match the size of data.
%
%       selection:: [] (required)
%           Data selection, defines which pixels should be taken into
%           account. The selection can be an empty matrix. In this case, it
%           is assumed that all pixels should be selected. The dimensions
%           of mask have to match the size of data.
%
% Example:
%   stxm_result = stxm(a_2d_diffraction_pattern,mask,[]);
%
% Output arguments:
%   result:: The integrated intensity of data.
%
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

    if isempty(sel)
        sel = ones(size(dat));
    end
    if isempty(mask)
        mask = zeros(size(dat));
    end

    result = sum(sum(data(~mask & sel)));
end