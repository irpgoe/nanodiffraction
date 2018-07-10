function [loc_r,loc_c,data_range] = find_xy_in_range(data,varargin)
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

    if nargin == 2        
        if isa(varargin{1},'function_handle')
            fh = varargin{1};
            data_range = fh(data);
        else
            lim = varargin{1};
            data_range = (data > lim(1)) & (data < lim(2));
        end
    else
        error('Please provide a lower and upper limit or a custom function handle');
    end
    
    [loc_r,loc_c] = ind2sub(size(data),find(data_range));
end