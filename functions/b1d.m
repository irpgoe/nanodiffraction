function [results] = b1d(dat,mask,sel,grid,varargin)
% B1D  performs an angular average based on histogramming the intensity of
% each pixel in a certain radial interval. B1D stands for *binning 1D* and
% is considered to be a rather simple but fast method to get accurate
% angular averages.
%   
%   result = b1d(data,mask,selection,grid,bins) 
%
% The following arguments are supported:
%     data:: [] (required)
%       The data that will be processed.
%
%     mask:: [] (required)
%       A logical mask that identifies all values that should not be taken 
%       into account during rebinning.
%
%     selection:: [] (required)
%       A logical array that selects all pixels that should be taken into
%       account. A pixel is taken into account if it is identified in the
%       selection as 1.
%
%     grid:: [] (required)
%       Usually, Qr is expected to be used, however, any radial grid can be
%       used here.
%
%     bins:: [360] (optional)
%       Number of bins.
%
% Example:
%   testavg = b1d(2d_diffraction_pattern,detector_mask,[],qr_map,512);
%
% Output arguments:
%   result:: Structure than contains the following fields:
%
%       - dat_1d:: One-dimensinoal azimuthally integrated data.
%
%       - qr:: X-axis, dependent on the grid used, however, it is usually
%       expected in units of the reciprocal wavevector qr (inv. nanometers)
%       and hence named qr.
%
%       - y:: A short-hand for dat_1d.
%
%    	- x:: A short-hand for qr. Also, independent of a specified grid.
%
%       - error:: Propagated measurement error.
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
	if nargin == 5
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
    
    % sort grid
    x = grid(~(mask | ~sel)); 
    [x, sortIndex] = sort(x);
    
    % binning of grid
    binEdge = linspace(min(x),max(x),bins);
    [n,~,bin] = histcounts(x,binEdge);
    binEdge = binEdge';
    
    % remove trailing zeros
    x = binEdge(1:sum(n~=0));
    
    % put intensity values into bins (each qr value has a bin associated) and take mean
    aavg = dat(~(mask | ~sel));
    aavg = aavg(sortIndex);
    aavg = accumarray(bin,aavg,[],@mean,NaN);
    
    % remove NaNs from intensity and bin index (NaNs appear where n is 0)
    n(isnan(aavg)) = [];
    aavg(isnan(aavg)) = [];
    
    % save results
    results.dat_1d = aavg;
    results.qr = x;
    results.y = aavg;
    results.x = x;
    results.error = sqrt(aavg)./sqrt(n');
end