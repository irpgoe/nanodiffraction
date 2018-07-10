function [result] = crystals(data,mask,sel,threshold)    
% CRYSTAL  Calculates Signal-to-Noise Ratio and
% summarizes reflection intensity (in between qr_min and qr_max if
% requested)
%
%   [RESULT] = GET_INDEX(DATA,MASK,SEL,THRESHOLD)
%
%   The following arguments are required:
%       DATA:: ()
%           Data (2d array).
%
%       MASK:: ()
%           Bad pixels.
%
%       SEL:: ()
%           Good pixels.
%
%       THRESHOLD:: ()
%           Maximal qr value to be taken into account.
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
        sel = ones(size(data));
    end
    if isempty(mask)
        mask = zeros(size(data));
    end
    
    data = data.*(~mask & sel);
    
    % several variations are possible, either
    result = sum(sum(data(data > threshold)));
    % or:
%     result = sum(sum(data((snr./sqrt(bgr)) > threshold)));
    % or:
%     result = sum(sum((snr./sqrt(bgr)) > threshold));
    % or:
%     result = sum(sum(snr((snr./sqrt(bgr)) > threshold)));
    
end