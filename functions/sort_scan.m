function out = sort_scan(data,rows,cols,varargin)
% SORT_SCAN  Description.
%   
%   result = sort_scan(data,rows,cols,scanmode)
%
% The following arguments are supported:
%     data:: [] (required)
%       One-dimensional or two-dimensional data array. If the data is
%       two-dimensional, each column will be treated seperately and stored
%       in the third dimension of the output array, since the first
%       dimension of the input array will be expanded into the first two
%       dimensions of the output array. Data can also be a cell array.
%
%     rows:: [] (required)
%       Number of scan rows.
%
%     cols:: [] (required)
%       Number of scan columns.
%
%     scanmode:: ['horz'] (required)
%            Several raster scanning modes are possible:
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
% Example:
%   Example missing.
%
% Output arguments:
%   result:: Description.
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