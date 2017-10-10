function [dpcx, dpcy] = CM_nano(Mass)
% R = [Spalte,Zeile]
% calculates center of mass for an image
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

    Mass = abs(Mass);

    M = sum(Mass(:));

    [I,J] = size(Mass);

    JJ = meshgrid(1:J,1:I);
    II = meshgrid(1:I,1:J).';

    R = [sum(sum(Mass.*JJ)); sum(sum(Mass.*II))];
    R = R/M;

    if I == 1
        R(2,1) = 0;
    end

    if J == 1
        R(1,1) = 0;
    end

    dpcx = R(1);
    dpcy = R(2);

end
