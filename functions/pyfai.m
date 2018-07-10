function result = pyfai(dat,mask,pby,pbz,varargin)
% PYFAI  performs an angular average based on the pyfai integration routine
% written by Jerome Kieffer, see 
% The fast azimuthal integration Python library: pyFAI,
% J. Kieffer, G. Ashiotis, A. Deschildre, Z. Nawaz, J. P. Wright, 
% D. Karkoulis, F. E. Picca 
% Journal of Applied Crystallography (2015) 48 (2), 510-519
%   
%   [RESULT] = PYFAI(DATA,MASK,PBY,PBZ,SETTINGS) 
%
%   The following options are required:
%
%     DATA::
%       The data that will be processed.
%
%     MASK::
%       A logical mask that sets all values that should not be taken into
%       account to NaN
%
%     PBY::
%       Primary beam positon (horizontal)
%
%     PBZ::
%       Primary beam positon (vertical)
%
%   The following options are optional:
%
%     SETTINGS:: 
%       Struct that can accept the following arguments
%           - angles: ()
%               ...
%           - phi1: ()
%               ...
%           - phi2: ()
%               ...
%           - cake: ()
%               ...
%           - to2D: ()
%               ...
%           - roi2D: ()
%               ...
%           - bins: ()
%               ...
%           - ai: ()
%               ...
%           - kwa_model: ()
%               ...
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

    if nargin < 4
        return;
    end
    if nargin >= 4

        % setup most simple detector with pixel unit length 1
        kwa_detector = pyargs(...
        'detector','Eiger4M',...
        'pixel1',1,...
        'pixel2',1);
        dim = py.tuple({size(dat,2),size(dat,1)});
        
        % set azimuthal integrator object
        ai = py.pyFAI.AzimuthalIntegrator(kwa_detector);
        ai.setFit2D(1,pbz,pby);

        % simple 1d transform without mask
        kwa_model = pyargs('correctSolidAngle',0,...
                           'method','numpy',...
                           'unit','r_mm',...
                           'error_model','poisson');
    
        % number of bins
        bins = 512;
        
        % mask
        if ~isempty(mask)
            reshaped_mask = py.numpy.reshape(reshape(mask,1,[]),...
                py.tuple({size(dat,2),size(dat,1)}));
            kwa_model = pyargs(...
                'correctSolidAngle',0,...
                'method','numpy',...
                'unit','r_mm',...
                'error_model','poisson',...
                'mask',reshaped_mask);
        end
        
        % default settings
        settings.cake = 'off';
        settings.to2D = 'off';
    end    
    if nargin == 5
        % overwrite default settings
        fields=fieldnames(varargin{1});        
        for i=1:size(fields,1)
            thisfield = fields{i};		
            settings.(thisfield) = varargin{1}.(thisfield);
        end

        if isfield(settings,'ai')
            ai = settings.ai;
        elseif isfield(settings,'kwa_model')
            kwa_model = settings.kwa_model;
        elseif isfield(settings,'bins')
            bins = settings.bins;
        end
    end 

    % reshape to fit
    dat = py.numpy.reshape(reshape(dat,1,[]),dim);

    % Dynamic cake Integration
    if strcmp(settings.cake,'on')
        if ~isempty(settings.angles)
            phi1 = settings.angles(akk,bkk) + settings.phi1;
            phi2 = settings.angles(akk,bkk) + settings.phi2;
            cake = py.tuple({phi1,phi2});
            kwa_model = pyargs('correctSolidAngle',0,...
                'method','numpy',...
                'unit','r_mm',...
                'error_model','poisson',...
                'azimuth_range',cake);
        end
    end

    if strcmp(settings.to2D,'off')
        % simple 1d mapping
        result_pyfai = ai.integrate1d(dat, bins, kwa_model);
        
        result.error = cell2mat(cell(result_pyfai{3}.tolist()));
        result.dat_1d = cell2mat(cell(result_pyfai{2}.tolist()));
        result.r = cell2mat(cell(result_pyfai{1}.tolist()))./(1*1e3);
    end
    if strcmp(settings.to2D,'on')
        % 2d mapping
        tmp = ai.integrate2d(dat, bins, kwa_model);
        tmp2 = py.numpy.reshape(tmp{1},py.tuple({360*bins,}));
        
        result.error = cell2mat(cell(tmp{3}.tolist()));
        result.dat_2d = reshape(cell2mat(cell(tmp2.tolist())),bins,360);
        
        if ~isempty(settings.roi2D)
            normalize_1d = sum(result.dat_2d(settings.roi2D(1):settings.roi2D(2),:)==0);
            result.dat_1d = sum(result.dat_2d(settings.roi2D(1):settings.roi2D(2),:)).*(settings.roi2D(2)-settings.roi2D(1)+1)./(settings.roi2D(2)-settings.roi2D(1)+1-normalize_1d);
%             result.dat_1d = sum(result.dat_2d(settings.roi2D(1):settings.roi2D(2),:));
        end
        result.r = cell2mat(cell(tmp{2}.tolist()))./(1*1e3); 
    end
end