
function [Ln] = symmetry(data, bgr, sigma, noise_level, obj)                   
% SYMMETRY  identifies peaks in the image and determines whether n-fold
% rotational symmetry exists.
%
%   LN = SYMMETRY(DAT,BGR,SIGMA,NOISELEVEL,OBJ)
%
%   The following arguments are required:
%       DAT:: ()
%           Diffraction pattern.
%
%       BGR:: ()
%           Background, that will be subtracted from diffraction pattern.
%
%       SIGMA:: ()
%           Parameter of Gaussian filter.
%
%       NOISELEVEL:: ()
%           SNR level that functions as a threshold for peak detection.
%
%       OBJ:: ()
%           Object that contains detector coordinates.
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

% hard-coded, theta window size
delta_theta = 20;

% signal-to-noise ratio
snr = abs(data-bgr)./sqrt(bgr);

% remove NaNs and Infs from snr
snr(isnan(snr)) = 0;
snr(isinf(snr)) = 0;

% Gaussian filter (remove noise)
h = fspecial('gaussian', [10 10], sigma);
h = h./sum(h(:));
snr = conv2(snr,h,'same');

% binary image (snr threshold)
binary_image = (snr >= noise_level);

% blob detection
CC = bwconncomp(binary_image);

% minimum number of objects for detection
if CC.NumObjects <= 1
    Ln = 0;
    return;
end

% calculate the centroids
S = regionprops(CC,'Centroid'); % output format (y,x)

% get orientation, reference axis is x axis and rotation is clockwise
for ii = 1:numel(S)
S(ii).pos = [obj.p.zaxis(round(S(ii).Centroid(1))) obj.p.yaxis(round(S(ii).Centroid(2)))];
phi = atan2(obj.p.zaxis(round(S(ii).Centroid(1))),obj.p.yaxis(round(S(ii).Centroid(2))));
S(ii).Orientation = phi*180/pi;
S(ii).Scale = sqrt((S(ii).pos(1))^2 + (S(ii).pos(2))^2);
S(ii).Size = 1;
end


% make unique pairs
pairs = cell(nchoosek(numel(S),2),1);
k = 0;
for i = 1:numel(S)
    remaining = 1:numel(S);
    remaining(i) = [];
    for j = remaining
        if j<i
            continue
        else
            k = k + 1; % index
            pairs{k} = struct('i',i,'j',j,'scale_i',S(i).Scale,'scale_j',S(j).Scale,'size_i',S(i).Size,'size_j',S(j).Size,'phi_i',S(i).Orientation,'phi_j',S(j).Orientation,'pos_i',S(i).pos,'pos_j',S(j).pos);
        end
    end
end

% determine the angle of rotation between pair
theta = zeros(numel(pairs));
distance = zeros(numel(pairs));
for ii = 1:numel(pairs)

    pair = pairs{ii};

    % (debugging) fprintf(1,'(%d,%d), theta: %f\n,',pair.i,pair.j,abs(pair.phi_i  - pair.phi_j));

    % if vectors are parallel, no rotational symmetry exists
    if abs(pair.phi_i  - pair.phi_j) < 10
        continue
    end

    % if radial distance between pair is too large,
    % skip
    if abs(pair.scale_i - pair.scale_j) > 10
        continue
    end

    % determine the angle of rotation between pair
    if abs(pair.phi_i - pair.phi_j) > 180
        theta(ii) = 360-abs(pair.phi_i - pair.phi_j);
    else
        theta(ii) = abs(pair.phi_i - pair.phi_j);
    end

    % determine distance of pair
    distance(ii) = (pair.scale_i + pair.scale_j)/2;

end
distance(theta == 0) = [];
theta(theta == 0) = [];

% calculate histogram of angles and rotation order
% distribution function
binEdges = 0:4:ceil(max(obj.p.R(:)));
[dist_histogram,~,bins] = histcounts(distance,binEdges);

Ln = zeros(1,length(dist_histogram)+1);
for ii = find(dist_histogram)

    % rotation histogram
    [rot_histogram,~] = histc(theta(bins == ii),1:180);
    if length(rot_histogram) ~= 180
        continue
    end
    response = order_estimation_function(rot_histogram,delta_theta);

    [~,Ln(ii+1)] = max(response);

end
response = order_estimation_function(histc(theta,1:180),delta_theta);
[~,Ln(1)] = max(response);
                        
end


function response = order_estimation_function(histogram,delta_theta)

    % angular response function
    response = zeros(11,1);
    for order = 2:12
        % hat function (window)
        g = ones(1,delta_theta);

        % rotations per order
        rot = 360/order;
        pts = rot:rot:180;
        x = zeros(1,180);
        x(round(pts)) = 1;

        % convolve hat function with rotations per order
        x = conv2(x,g,'same');

        % the same for out-of-phase rotation
        offset_pts = rot/2:rot:180;
        x_offset = zeros(1,180);
        x_offset(round(offset_pts)) = 1;
        x_offset = conv2(x_offset,g,'same');

        % calculate response
        response(order) = (sum(x.*histogram) - sum(x_offset.*histogram))/order;
    end
    response(response<0) = 0;
    response = response/max(response(:));
end
