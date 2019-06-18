function [idx_out,closest] = get_index(x,val)
% GET_INDEX  calculates the index of a data entry within the data array which is closest to a reference value
%
%   ``[IDX,Y] = GET_INDEX(X,VAL)``
%
% Parameters
% ----------
% X: One-dimensional numerical array
%   Data array
%
% VAL: Integer of floating-point value or array of integer or floating
%   point values
%   Reference value or values
%
% Returns
% -------
% IDX: Integer
%   Index within the search array Y
% 
% Y: Integer of floating-point value
%   Value of the entry at position I in the search array Y
%
% Notes
% -----
% Example 1:
%
% .. code-block:: matlab
%
%   x = 1:100;
%   [idx,y] = get_index(x,22.7);
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
    
    idx_out = zeros(numel(val),1);
    closest = zeros(numel(val),1);
    for i = 1:numel(val)
        [~,idx_out(i)] = min(abs(x-val(i)));
        closest(i) = x(idx_out(i));
    end
end