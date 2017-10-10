function [results] = b1d(dat,mask,sel,grid,varargin)
% B1D  performs an angular average based on histogramming the intensity of
% each pixel in a certain radial interval. B1D stands for *binning 1D* and
% is considered to be a rather simple but fast method to get accurate
% angular averages.
%   
%   [RESULT] = B1D(DATA,sel,GRID,BINS) 
%
%   The following options are required:
%
%     DATA::
%       The data that will be processed.
%
%     MASK::
%       A logical sel that sets all values that should not be taken into
%       account to NaN
%
%     SEL::
%       A logical sel that selects all pixels that should be taken into
%       account 
%
%     GRID::
%       Usually, Qr is expected to be used, however, any radial grid can be
%       used here.
%
%   The following options are optional:
%
%     BINS:: (360)
%       Number of bins.
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
            
    % Settings, if none are passed than default values are taken
	if (nargin == 5)
		bins = varargin{1};
    else
        bins = 360;
    end
    if isempty(sel)
        sel = ones(size(dat));
    end
    if isempty(mask)
        mask = zeros(size(dat));
    end
    
    % apply sel on radius in pixels and data
    grid(mask | ~sel) = NaN;
    dat(mask | ~sel) = NaN;

    % transform into one-dimensional array
    aavg = squeeze(reshape(dat,1,[]));    
    x = squeeze(reshape(grid,1,[]));

    % remove NaNs
    aavg(isnan(aavg)) = [];
    x(isnan(x)) = [];
    
    % sorting (x_i, I_i) pairs
    [x, sortIndex] = sort(x);
    aavg = aavg(sortIndex);
    
    % binning of x values
    binEdge = linspace(min(x),max(x),bins);
    [n,~,bin] = histcounts(x,binEdge);
    
    % put intensity values into bins (each qr value has a bin associated) and take mean
    aavg = accumarray(bin',aavg,[],@mean,NaN);
    binEdge = binEdge';
    
    % remove NaNs from intensity
    aavg(isnan(aavg)) = [];

    % binning and remove trailing zeros
    x = binEdge(1:sum(n~=0));
    
    % save results
    results.dat_1d = aavg;
    results.qr = x;
    results.error = sqrt(aavg);
end