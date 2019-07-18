function out = sort_scan(data,rows,cols,varargin)
% SORT_SCAN  Reshapes a given one-dimensional or two-dimensional numerical array or cell array into the dimensions of the scan. In the case of a two-dimensional array, each column is treated independently and the result will be output as a stack of two-dimensional numerical arrays. 
%   
%   ``out_array = sort_scan(data,rows,cols,scanmode)``
%
% Parameters
% ----------
%   data: One-dimensional numerical array or cell array
%       One-dimensional or two-dimensional numerical array. If the data is
%       two-dimensional, each column will be treated seperately and stored
%       in the third dimension of the output array, since the first
%       dimension of the input array will be expanded into the first two
%       dimensions of the output array. Input data can also be in form of a cell array.
%
%   rows: Integer
%       Number of scan lines
%
%   cols: Integer
%       Number of scan points per line
%
%   scanmode: Default = 'horz'
%       Several raster scanning modes are possible, see ASCII art in the following code block:
%
%       .. code-block:: matlab
%
%             1 horizontal <horz>
%                *===================>
%                ====================>
%                ========>...            
%             2 horizontal, alternating <horzalt>
%                *===================>
%                                    v
%                <===================<
%                v
%             3 vertical <vert>
%                *  v  v  v  v
%                v  v  v  v  .
%                v  v  v  v  .
%                v  v  v  v  .
%             4 vertical,alternating <vertalt>
%                *   A > v   A > v
%                v   A   v   A   .
%                v   A   v   A   .
%                v > A   v > A   .
%
% Returns
% -------
%   out_array: Two- or three-dimensional numeric array
%       Two- or three-dimensional numeric array. The first two dimensions are of the size of the scan dimensions.
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

    type = 'horz';
    if nargin > 3
        type = varargin{1};
    end

    % data input type
    if iscell(data)
        out = cell(rows,cols,size(data,2));
    elseif isstruct(data)
        out = repmat(data(1),rows,cols,size(data,2));
    else
        out = zeros(rows,cols,size(data,2));
    end

    % transform row to column
    if size(data,1) == 1
        data = data';
    end
    
    % treat columns seperately
    for ii = 1:size(data,2)
        data_per_col = data(:,ii);
        switch type
            case 'horz'
                data_per_col = transpose(reshape(data_per_col,[cols rows]));
            case 'horzalt'
                data_per_col = transpose(reshape(data_per_col,[cols rows]));
                data_per_col(2:2:end,:) = fliplr(data_per_col(2:2:end,:));
            case 'vert'
                data_per_col = reshape(data_per_col,[rows cols]);
            case 'vertalt'
                data_per_col = reshape(data_per_col,[rows cols]);
                data_per_col(:,2:2:end) = flipud(data_per_col(:,2:2:end));
            otherwise
                error('Choose between <horz>,<horzalt>,<vert>,<vertalt>');
        end
        out(:,:,ii) = data_per_col;
    end
end