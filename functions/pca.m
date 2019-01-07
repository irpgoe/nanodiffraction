function [results] = pca(dat,qy,qz,mask,sel)
% PCA  calculates the covariance of the diffraction pattern and determines the eigenvalues of the scattering distribution.
%
%   result = pca(data,qy,qz,mask,selection)
%
% The following arguments are supported:
%       data:: [] (required)
%           Diffraction pattern.
%
%       qy:: [] (required)
%           2d array of horizontal momentum transfer. The array has to
%           match the dimensions of data.
%
%       qz:: [] (required)
%           2d array of vertical momentum transfer. The array has to match 
%           the dimensions of data.
%
%       mask:: [] (required)
%           Mask that defines which pixels are considered bad. The mask can
%           be an empty matrix. In this case, it is assumed that all pixels
%           are valid. The array has to math the dimensions of data.
%
%       selection:: [] (required)
%           Data selection matrix. The selection can be an empty matrix. In
%           this case, it is assumed that all pixels should be selected. 
%           The array has to match the dimensions of data.
%
% Example:
%   pca_result = pca(a_2d_diffraction_pattern,qy,qz,mask,[]);
%
% Output arguments:
%   result:: A structure that contains the following fields:
%       - w:: Anisotropy of the data.
%
%       - angle:: Orientation angle of the largest component.
%
%       - lambda_1:: Eigenvalue of the largest component.
%
%       - lambda_2:: Eigenvalue of the smallest component.
%
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

    % apply (pca) mask
    if isempty(mask)
        mask = zeros(size(dat));
    end
    if isempty(sel)
        mask = ones(size(dat));
    end
    dat = abs(dat.*(~mask & sel));
    
    % first order moments
    tr = sum(dat(:));
    
    % calculate covariance
    cov_qy_qy=sum(sum(dat.*(qy.^2)))./tr;
    cov_mix=sum(sum(dat.*qy.*qz))./tr; 
    cov_qz_qz=sum(sum(dat.*(qz.^2)))./tr;

    % covariance matrix
    A=[cov_qy_qy cov_mix; cov_mix cov_qz_qz];

    % solve eigenvalue problem (V: Eigenvectors, D:
    % diagonal matrix)
    % Eig must not contain NaN or Inf, hence
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    [V,D] = eig(A);

    % Orientation of largest Eigenvector
    if D(1,1)>=D(2,2)
        angle = atan2d(V(2,1),V(1,1));
        w = (D(1,1)-D(2,2))/(D(1,1)+D(2,2));
        lambda_1 = D(1,1);
        lambda_2 = D(2,2);
    else
        angle = atan2d(V(2,2),V(1,2));
        w = (D(2,2)-D(1,1))/(D(1,1)+D(2,2));
        lambda_1 = D(2,2);
        lambda_2 = D(1,1);
    end    
    
    % Project onto range [0 180)
    if angle < 0
        angle = angle + 180;
    end
    
    % use mex function
%     if ~isempty(empty)
%     [w,angle] = mex_pca(dat-empty,mask,qy,qz);        
%     else
%     [w,angle] = mex_pca(dat,mask,qy,qz);
%     end

    
    % Save results
    results.w = w;
    results.angle = angle;
    results.lambda_1 = lambda_1;
    results.lambda_2 = lambda_2;
end