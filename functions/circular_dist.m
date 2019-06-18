function [circMean, circVar] = circular_dist(phi,I)
% CIRCULAR_DIST  Calculates the circular mean and variance of a 
% distribution based on 
% https://en.wikipedia.org/wiki/Directional_statistics
%
%   ``[circMean, circVar] = circular_dist(phi,I)``
%
% In a principal component analysis of a scattering distribution, the covariance matrix of the
% distribution is calculated and orthogonalized. The procedure is described in more detail in
% Bernhardt et al., 2016.
%
% Parameters
% ----------
% phi: One-dimensional numeric array
%   phi contains angles within the range [-180, 180]
%
% I: One-dimensional numeric array
%   Intensity as a function of the azimuthal angle phi
%
% Returns
% -------
% circMean: Numeric value
%   The circular mean of the distribution
%
% circVar: Numeric value
%   The circular variance of the distribution
%
% Notes
% -----
% Example 1:
%
% .. code-block:: matlab
%
%       e = nanodiffraction();
%       streak = e.simulate_streak();
%       [phi,I] = split_struct(b1d(streak,[],[],e.phi,360),{'x','y'});
%       [circMean, circVar] = circular_dist(phi,I);
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

    % correct phi
    phi = phi*2-180;
    
    % number of data points
    N = numel(phi);
    
    % transform to complex number
    z = zeros(N,1);
    for jj = 1:N
        z(jj) = I(jj)*exp(1i*phi(jj)*pi/180);
    end

    % circular mean
    circMean = angle(sum(z)/N)*180/pi;
    
    % correct circ. mean
    circMean = (circMean+180)/2;
    
    % circular variance
    c = sum(cos(angle(z)));
    s = sum(sin(angle(z)));
    r = sqrt(c^2 + s^2)/N;
    circVar = 1-r; 
end    