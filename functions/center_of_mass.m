function [dpcy, dpcz] = center_of_mass(data)
% CENTER_OF_MASS calculates the center of mass of a diffraction pattern.
% This is commonly used for experiments without beamstop where only the
% shift of the center of mass of the primary beam due to refraction of the
% beam is recorded. When applied in a raster-scan experiment, one obtains a
% differential phase contrast image of the sample.
%
%   ``[dpcy, dpcz] = center_of_mass(data)``
%
% Parameters
% ----------
% data: numerical array
%   Diffraction pattern
%
% Returns
% -------
% dpcy: Numeric value
%   center of mass along the y-axis
%
% dpcz: Numeric value
%   center of mass along the z-axis
%
% Notes
% -----
% Example 1:
%
% .. code-block:: matlab
%
%       example = missing;
%
% See also 

% Copyright 2017 Institute for X-ray Physics (University of Göttingen)
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
% associated documentation files (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge, publish, distribute,
% sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or
% substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
% NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
% OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    data = abs(data);

    M = sum(data(:));

    [I,J] = size(data);

    JJ = meshgrid(1:J,1:I);
    II = meshgrid(1:I,1:J).';

    R = [sum(sum(data.*JJ)); sum(sum(data.*II))];
    R = R/M;

    if I == 1
        R(2,1) = 0;
    end

    if J == 1
        R(1,1) = 0;
    end

    dpcy = R(1);
    dpcz = R(2);

end
