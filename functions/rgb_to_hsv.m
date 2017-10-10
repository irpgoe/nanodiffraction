function [hue,sat,val,cmap] = rgb_to_hsv(rgb)
% RGB_TO_HSV  Transform RGB values to HSV values. 
%
%   [HUE,SAT,VAL,CMAP] = RGB_TO_HSV(RGB)
%
%   The following arguments are required:
%       RGB:: ()
%           Array with three entries [R G B].
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

    % seperate entries
    R = rgb(1);
    G = rgb(2);
    B = rgb(3);
    
    [Vmax,Imax] = max([R G B]);
    [Vmin,Imin] = min([R G B]);

    % siehe hierzu https://de.wikipedia.org/wiki/HSV-Farbraum

    % hue
    if Vmax == Vmin
        H = 0;
    elseif Imax == 1
        H = 60*(0 + (G-B)/(Vmax-Vmin));
    elseif Imax == 2
        H = 60*(2 + (B-R)/(Vmax-Vmin));
    elseif Imax == 3
        H = 60*(4 + (R-G)/(Vmax-Vmin));
    end

    if H < 0
        H = H + 360;
    end

    % sat
    if Vmax == 0
        S = 0;
    else
        S = (Vmax - Vmin) / Vmax;
    end

    % val
    V = Vmax;

    hue = H;
    sat = S;
    val = V;

    cmap(1) = hue;
    cmap(2) = sat;
    cmap(3) = val;
end

    