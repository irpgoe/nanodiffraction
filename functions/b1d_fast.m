function [results] = b1d_fast(dat,mask,sel,sortIndex,n,bin)
% B1D  performs an angular average based on histogramming the intensity of
% each pixel in a certain radial interval. B1D stands for *binning 1D* and
% is considered to be a rather simple but fast method to get accurate
% angular averages.
%   
%   [RESULT] = B1D(DATA,MASK,SEL,GRID,BINS) 
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
    if isempty(sel)
        sel = ones(size(dat));
    end
    if isempty(mask)
        mask = zeros(size(dat));
    end
    
    % put intensity values into bins (each qr value has a bin associated) and take mean
    aavg = dat(~(mask | ~sel));
    aavg = aavg(sortIndex);
    aavg = accumarray(bin,aavg,[],@mean,NaN);
    
    % remove NaNs from intensity and bin index (NaNs appear where n is 0)
    n(isnan(aavg)) = [];
    aavg(isnan(aavg)) = [];
    
    % save results
    results.dat_1d = aavg;
    results.error = sqrt(abs(aavg))./sqrt(n');
end