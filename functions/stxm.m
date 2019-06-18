function result = stxm(data,mask,sel)
% STXM  Integrated the intensity in a diffraction pattern or part of a diffraction pattern
%
%   ``result = stxm(data,mask,selection)``
%
% Parameters
% ----------
%   data: two-dimensional numeric array 
%       A diffraction pattern 
%
%   mask: logical array, default = []
%       Logical array of the size of the input data that defines which pixels are considered bad (1 = bad, 0 = valid)
%
%   selection: logical array, default = []
%       Logical array of the size of the input data that defines which pixels should be analyzed (1 = valid, 0 = invalid)
%
% Returns
% -------
%   result: Double 
%       The integrated intensity of data within a given selection and without masked pixels
%
% Notes
% -----
% Example 1:
%   
% .. code-block:: matlab
%
%   stxm_result = stxm(data,[],[]);
%
% See also PCA

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
        sel = ones(size(data));
    end
    if isempty(mask)
        mask = zeros(size(data));
    end

    result = sum(sum(data(~mask & sel)));
end