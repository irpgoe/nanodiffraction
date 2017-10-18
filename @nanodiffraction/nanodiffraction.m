classdef nanodiffraction < handle
    % NANODIFFRACTION  Class that defines the experimental geometry and
    % calculates geometrical parameters accordingly. 
    % The key method is "analyze_scan" that is used to process x-ray 
    % diffraction data obtained from scanning experiments.
    %
    %   EXPERIMENT = NANODIFFRACTION(ARGLIST) creates an object, an
    %   experiment, with a list of key-value pairs provided in ARGLIST.
    %   ARGLIST consists of NAME - VALUE pairs. NAME is a string describing 
    %   the option and VALUE can be any input. See help inputParser for 
    %   more information and an example below. Typically, a specific input 
    %   is requested by default and an error is thrown if the types do not
    %   agree.
    %
    %   The nanodiffraction class holds two types of parameters,
    %   experiment-specific parameters and scan parameters. In
    %   addition, it contains smaller useful helper functions such as 
    %   q_of_n() that takes in pixel coordinates and transforms them into q 
    %   coordinates.
    %   In order to perform calculations on real data, one has to couple
    %   the nanodiffraction class to the files class to "tell the
    %   experiment, where the data is stored". This can be achieved by
    %   setting EXPERIMENT.data = [handle to files class] or by using the
    %   helper function link([handle to files class],[handle to helper
    %   class]). See help link for more information.
    %
    %   The following options are supported:
    %
    %   energy:: [0] (not necessary if wavenumber or wavelength are given)
    %     x-ray photon energy in [eV]
    %
    %   wavenumber:: [0] (not necessary if energy or wavelength are given)
    %     x-ray photon wavenumber in [m^-1]
    %
    %   wavelength:: [0] (not necessary if wavenumber or energy are given)
    %     x-ray photon wavelength in [m]
    %
    %   detDistance:: [5]
    %     detector distance in [m]
    %
    %   pixelsize:: [172E-6]
    %     detector pixel size in [m]
    %
    %   Ny:: [487]
    %     number of detector pixels along y axis (horizontal axis)
    %
    %   Nz:: [619]
    %     number of detector pixels along z axis (vertical axis)
    %
    %   pby:: [1]
    %     horizontal primary beam position relative to first pixel in 
    %     top-left corner (when viewed along beam propagation direction)
    %
    %   pbz:: [1]
    %     vertical primary beam position relative to first pixel in 
    %     top-left corner (when viewed along beam propagation direction)
    %
    %   fluoCounters:: [{}]
    %     Cell array of thresholds to generate maps of different element 
    %     distributions    
    %
    %   data:: [{}]
    %     files-class handle.
    %
    %   Example::
    %      experiment =
    %      nanodiffraction('energy',12.8e3,'detDistance',1.2,...
    %      'pby',121.2,'pby',144.1,...
    %      'Ny',1100,'Nz',1220,...
    %      'fluoCounters',{{'total',1 4100},{'caka',600,660}});
    %
    %   Properties::
    %      <experiment>  
    %                   - scan (holds scan parameters)
    %                   - data (links to the data)
    %                   - helper (stores helper functions)
    %
    %   Methods::
    %      See seperate help for each method.
    
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

    properties
        helper = struct;    % holds helper functions
        scan = struct;      % holds scan parameters
        data = files;       % holds link to the data
        
        % experimental parameters
        energy = 0;
        wavelength = 0;
        wavenumber = 0;
        yaxis = [];
        zaxis = [];
        R = [];
        theta = [];
        Ny = 0;
        Nz = 0;
        pby = 0;
        pbz = 0;
        pixelsize = 0;
        detDistance = 0;
        qyaxis = [];
        qzaxis = [];
        qr = [];
        phi = [];
        
        % wait parameters used for data polling
        wait = 0;
        wait_it = 0;
        wait_maxit = 1000;
        
        % original values that are overwritten once a detector roi or a
        % binning are chosen
        Ny_orig = 0;
        Nz_orig = 0;
        pby_orig = 0;
        pbz_orig = 0;
        pixelsize_orig = 0;
        
        % mask, selections and corrections
        mask = [];
        sels = [];
        corr = [];
        mask_orig = [];
        sels_orig = [];
        corr_orig = [];
        
        % empty
        empty = [];
        
        % roi and binning
        detectorRoi = 'off';
        binning = 'off';
        binx = 1;
        biny = 1;
        detRoiY = [];
        detRoiZ = [];
        
        % fluorescence counters
        fluoCounters = cell(1);
    end
    methods
        % constructor
        function obj = nanodiffraction(varargin)
            
        	pin = inputParser;
            % set defaults and add optional arguments
            addOptional(pin,'energy',13.8);
            addOptional(pin,'wavenumber',0);
            addOptional(pin,'wavelength',0);
            addOptional(pin,'detDistance',5);
            addOptional(pin,'pixelsize',172E-6);
            addOptional(pin,'Ny',487);
            addOptional(pin,'Nz',619);
            addOptional(pin,'pby',1);
            addOptional(pin,'pbz',1);
            addOptional(pin,'data',files);
            addOptional(pin,'fluoCounters',{'total',60,2300},@iscell);
            
            % parse arguments
            parse(pin,varargin{:});
            names = fieldnames(pin.Results);
            
            % copy data to class property and all other parameters to
            % the parameter struct (p)
            for ii = 1:numel(names)
                if strcmp(names{ii},'data')
                    obj.data = pin.Results.data;
                else
                    obj.(names{ii}) = pin.Results.(names{ii});
                end
            end

            % calculate missing energy, wavelength or wavenumber
            if obj.energy==0 && obj.wavelength==0 && obj.wavenumber==0
                warning('Please specify an energy, wavenumber or wavelength.');
            else
                hbar = 0.19732697E-6;
                if obj.energy == 0
                    % energy is not set
                    if obj.wavelength == 0 && obj.wavenumber ~= 0
                        % wavenumber is set
                        obj.energy = hbar * obj.wavenumber;
                        obj.wavelength = 2*pi / obj.wavenumber;
                    else
                         % wavelength is set
                        obj.energy = hbar * 2*pi / obj.wavelength;
                        obj.wavenumber = 2*pi / obj.wavelength;
                    end

                else
                    % energy is set
                    if obj.wavelength == 0
                        % calculate wavelength from energy
                        obj.wavelength = hbar * 2*pi / obj.energy;
                    end
                    if obj.wavenumber == 0
                        % calculate wavenumber from energy  
                        obj.wavenumber = obj.energy / hbar;
                    end
                end
            end
            
            % helper functions
            % get qr from pixel number or pixel number from qr
            n_of_q = @(q) tan(2*asin(q*obj.wavelength*1E9/4/pi)) * obj.detDistance/obj.pixelsize;
            q_of_n = @(n) 4*pi/obj.wavelength/1E9 * sin(atan(n*obj.pixelsize/obj.detDistance)/2);
            % get two_theta from pixel number
            twoTh_of_n = @(n) atand(n.*obj.pixelsize./obj.detDistance);
            % transform a 2d array into a column vector
            to1dcol = @(data) reshape(data,[],1); % colwise -> col
            % transform a 2d array into a row vector
            to1drow = @(data) reshape(data',1,[]); % rowwise -> row
            % transform a 1d array into a 2d matrix with dimensions [ny nx]
            % (works only for column vector)
            to2d = @(data,ny,nz) reshape(data,[nz ny]); 
            % description missing
            reduce = @(i,start,skip) start + (i-1)*skip;
            % split up struct containing dat_1d and qr (from azimuthal 
            % averaging) into seperate variables
            simplify_sf = @(data) deal(1:length(data.dat_1d),data.dat_1d);
            % rescale distribution to [0 1)
            scale = @(amp) (amp - min(min(amp)))/(max(max(amp)) - min(min(amp)));
            scaleH = @(amp,high) scale(amp).*high;
            scaleLH = @(amp,low,high) scaleH(amp,high-low) + low;
            % clip distribution to [low high]
            clipL = @(amp,low) amp.*(1-(amp<low)) + low*(amp<low);
            clipH = @(amp,high) amp.*(1-(amp>high)) + high*(amp>high);
            clipLH = @(amp,low,high) clipL(clipH(amp,high),low);
            % calculate 5% and 95% percentile from distribution or map
            autoClip = @(amp) deal(prctile(reshape(amp,1,[]),5),prctile(reshape(amp,1,[]),95));
            % find index of element at position pos (in arb. units)
            getIndex = @(qr,pos) min(abs(qr-pos));
            % add transparency
            addTransp = @(transp) set(imhandles(gcf),'AlphaData',transp);
            % normalize distribution
            normDist = @(im) (im - mean(reshape(im,1,[])))/(std(reshape(im,1,[])));
            
            % save helper to struct
            obj.helper.n_of_q = n_of_q;
            obj.helper.q_of_n = q_of_n;
            obj.helper.twoTh_of_n = twoTh_of_n;
            obj.helper.to1dcol = to1dcol;
            obj.helper.to1drow = to1drow;
            obj.helper.to2d = to2d;
            obj.helper.reduce = reduce;
            obj.helper.simplify_sf = simplify_sf;
            obj.helper.scale = scale;
            obj.helper.scaleH = scaleH;
            obj.helper.scaleLH = scaleLH;
            obj.helper.clipL = clipL;
            obj.helper.clipH = clipH;
            obj.helper.clipLH = clipLH;
            obj.helper.autoClip = autoClip;
            obj.helper.getIndex = getIndex;
            obj.helper.addTransp = addTransp;
            obj.helper.normDist = normDist;
            
            % detector x-y-axis 
            obj.yaxis = (1:obj.Ny) - obj.pby;
            obj.zaxis = obj.pbz - (1:obj.Nz);
            
            % detector r-theta-axis
            [xx,yy] = meshgrid(obj.yaxis,obj.zaxis);
            obj.R = sqrt(xx.^2 + yy.^2);
            obj.theta = atand(obj.R.*obj.pixelsize./obj.detDistance);
            obj.phi = atan2d(repmat(obj.zaxis',1,obj.Ny),repmat(obj.yaxis,obj.Nz,1));
            obj.phi(obj.phi<0) = obj.phi(obj.phi<0) + 360;
            
            % detector q-axis
            obj.qyaxis = obj.helper.q_of_n(obj.yaxis);
            obj.qzaxis = obj.helper.q_of_n(obj.zaxis);
            obj.qr = obj.helper.q_of_n(obj.R);
            
            % initialize detector roi
            obj.detRoiY = [1 obj.Ny];
            obj.detRoiZ = [1 obj.Nz];
            
            % detector mask
            obj.mask = zeros(obj.Nz,obj.Ny);
            obj.mask_orig = zeros(obj.Nz,obj.Ny);
            
            % empty
            obj.empty = zeros(obj.Nz,obj.Ny);
            
            % detector roi and binning off
            obj.detectorRoi = 'off';
            obj.binning = 'off';
            obj.binx = 1;
            obj.biny = 1;
            
            % initial PCA mask
            obj.sels = ones(obj.Nz,obj.Ny);
            obj.sels_orig = ones(obj.Nz,obj.Ny);
            
            % initial correction
            obj.corr = ones(obj.Nz,obj.Ny);
            obj.corr_orig = ones(obj.Nz,obj.Ny);
            
            % save the original values for the detector before binning or
            % anything else that might affect number of pixels, primary
            % beam position or pixel size
            obj.Ny_orig = obj.Ny;
            obj.Nz_orig = obj.Nz;
            obj.pby_orig = obj.pby;
            obj.pbz_orig = obj.pbz;
            obj.pixelsize_orig = obj.pixelsize;
        end
        
        
        function [n] = n_of_q(obj,q)
            % N_OF_Q  Calculates pixel coordinate in pixels given a Q
            % value.
            %
            %   N = N_OF_Q(Q) where Q is the Q coordinate in [nm^-1]
            %
            % This function cannot be defined as an anonymous function,
            % because pixelsize might in fact not be constant (and the
            % other constants in some rare occasions as well)
            
            n = tan(2*asin(q*obj.wavelength*1E9/4/pi)) * obj.detDistance/obj.pixelsize;
        end
        
        function [q] = q_of_n(obj,n)
            % Q_OF_N  Calculates Q coordinate in [nm^-1] given a pixel
            % coordinate.
            %
            %   Q = Q_OF_N(N) where N is the pixel coordinate
            %
            % This function cannot be defined as an anonymous function,
            % because pixelsize might in fact not be constant (and the
            % other constants in some rare occasions as well)
            
            q = 4*pi/obj.wavelength/1E9 * sin(atan(n*obj.pixelsize/obj.detDistance)/2);
        end
        
        function refresh(obj)
            % REFRESH  This function might be needed in case the geometry 
            % settings change. This could happen eg. when for 
            % some reason the detector pixel size was changed (e.g. a roi 
            % was chosen) and q_r, primary beam position, etc. has to be 
            % recalculated.
            %
            %   REFRESH()
                
        
            fprintf(1,'Resetting detector grid...\n');
            
            % detector x-y-axis 
            obj.yaxis = (1:obj.Ny) - obj.pby;
            obj.zaxis = obj.pbz - (1:obj.Nz);
            
            % detector r-theta-axis
            [xx,yy] = meshgrid(obj.yaxis,obj.zaxis);
            obj.R = sqrt(xx.^2 + yy.^2);
            obj.theta = atand(obj.R.*obj.pixelsize./obj.detDistance);
            
            % detector q-axis
            obj.qyaxis = obj.q_of_n(obj.yaxis);
            obj.qzaxis = obj.q_of_n(obj.zaxis);
            obj.qr = obj.q_of_n(obj.R);
            obj.phi = atan2d(repmat(obj.zaxis',1,obj.Ny),repmat(obj.yaxis,obj.Nz,1));
            obj.phi(obj.phi<0) = obj.phi(obj.phi<0) + 360;
            
            % update masks
            if max(size(obj.mask) ~= [obj.Nz,obj.Ny])
                fprintf(1,'Processing detector mask...\n');
                obj.mask = obj.process_mask(obj.mask);
            end
            if max(size(obj.sels(:,:,1)) ~= [obj.Nz,obj.Ny])
                fprintf(1,'Processing selections...\n');
                tmp = zeros(obj.Nz,obj.Ny,size(obj.sels,3));
                for ii = 1:size(obj.sels,3)
                    tmp(:,:,ii) = obj.process_mask(obj.sels(:,:,ii));
                end
                obj.sels = tmp;
            end     
            if max(size(obj.corr) ~= [obj.Nz,obj.Ny])
                fprintf(1,'Processing correction matrix...\n');
                obj.corr = obj.process(obj.corr);
            end
            
            
        end
        
        
        function new = copy(obj)
            % COPY  Makes a copy of a handle object. This keeps track of 
            % all the given parameters so that multiple instances of an 
            % experiment can be used with only slight changes in 
            % experimental setting.
            %
            %   NEW_INSTANCE = COPY(NANODIFFRACTION_OBJECT)
            %
            %   EXAMPLE::
            %       waxs = copy(generic_experiment)
            %       waxs.p.detDistance = 0.2;
            %       saxs = copy(generic_experiment)
            %       saxs.p.detDistance = 5;
            
            
            % Instantiate new object of the same class.
%             fprintf(1,'/* (ignore) ');
            new = feval(class(obj));
%             fprintf(1,'\b */\n');
            
            % Copy all non-hidden properties.
            prop = properties(obj);
            for i = 1:length(prop)
                new.(prop{i}) = obj.(prop{i});
            end
        end
        
        
        function result = analyze_scan(obj,varargin)
            % ANALYZE_SCAN  Main method of nanodiffraction class. It
            % analyzes a scan based on several implemented methods that can
            % be chosen in a random order and combination.
            %
            %   RESULT = ANALYZE_SCAN(VARARGIN) output the results as a
            %   struct where each entry in the struct is named after the
            %   method that was chosen for analysis. All relevant data can
            %   be found therein.
            %
            %   The following arguments are accepted:
            %
            %      roi:: []
            %        A logical mask of size SNy x SNz. A logical (1)
            %        indicates a scan point, that will be analyzed.
            %
            %      yCrop, zCrop:: [start end] 
            %        start and end can be given to define a certain
            %        range within the scan to be analyzed
            % 
            %      ySkip, zSkip:: [1]
            %        If only every n-th scan point should be analyzed
            %
            %      emptySub:: [off] 
            %        Can be switched between 'on' | 'off' . Subtracts the
            %        empty image stored in <nanodiffraction>.scan.empty
            %        from the data immediately after the data has been
            %        read.
            %
            %      correction:: [off] 
            %        Can be switched between 'on' | 'off' . Multiplies a 2d
            %        correction matrix with the data. E.g. used to account
            %        for semitransparent objects obscuring the detector
            %
            %      method:: [] 
            %        A string composed of the following options: 'debug', 
            %        'rotate','pca','stxm','maxProj','crystal','fluo',
            %        'average','sum','data'. Methods can be combined using
            %        a '+' symbol between the methods:
            %           e.g. 'stxm+pca+maxProj+sum'
            %
            %      fluoCounters:: [{}] 
            %         The format is {{<name>,<channel_low>,<channel_high>},...}, 
            %           e.g. {{'caka',1000,1400},{},...}.
            %
            %      parallel:: [off]
            %        When turned on, all available workers will be used to
            %        speed up the calculation process. This is currently
            %        only implemented for the stxm method.
            %
            %      scanmode:: [horz|horzalt|vert|vertalt]
            %        See description below.
            %
            % Some methods can be given a set of arguments in a struct.
            % Note, that all struct parameters are optional. If some or all
            % arguments are not given, then default values will be used, as
            % indicated in paranthesis [default value].
            %
            %   pca: Performs a PCA analysis as e.g. demonstrated in 
            %   Bernhardt, et.al., X-Ray Micro- and Nanodiffraction Imaging 
            %   on Human Mesenchymal Stem Cells and Differentiated Cells,
            %   Biophys. J. (2016), 110, 680-690.
            %
            %       PARAMETER ARGUMENTS: 'pcaArgs',struct()
            %           struct() can take the following arguments:
            %           - 'qy': [repmat(obj.qyaxis ,obj.Nz,1)]
            %           - 'qz': [repmat(obj.qzaxis',1,obj.Ny)]
            %           - 'sel': [obj.sels]
            %           - 'mask': [obj.mask]
            %
            %   stxm: Integrates the scattered intensity of each
            %   diffraction pattern (darkfield contrast).
            %
            %       PARAMETER ARGUMENTS: 'stxmArgs',struct()
            %           struct() can take the following arguments:
            %           - 'mask':   [obj.mask]
            %           - 'qrMin':  [-1]
            %           - 'qrMax':  [-1]
            %           - 'sel':    [ones(obj.Nz,oabj.Ny))]
            %
            %   maxProj: Calculates the maximum projection of a stack of
            %   diffraction patterns. No arguments required.
            %
            %   crystal: Performs a thresholding for crystal mapping based
            %   on SNR criterium.
            %
            %       PARAMETER ARGUMENTS: 'crystalArgs',struct()
            %           struct() can take the following arguments:
            %           - Description missing
            %           -
            %           -
            %           -            
            %
            %   pyfai: Fast azimuthal integration based on the python pyfai
            %   library.
            %
            %       PARAMETER ARGUMENTS: 'pyfaiArgs',struct()
            %           struct() can take the following arguments:
            %           - Description missing
            %           -
            %           -
            %           -         
            %
            %   b1d: Fast azimuthal integration based on rebinning the
            %   data.
            %
            %       PARAMETER ARGUMENTS: 'b1dArgs',struct()
            %           struct() can take the following arguments:
            %           - 'grid':   [obj.qr]
            %           - 'bins':   [100]
            %           - 'mask':   [obj.mask]
            %           - 'sel':    [ones(obj.Nz,obj.Ny)]
            %
            %   symmetry: Identifies reflections in the data at calculates
            %   rotational order. Used e.g. in Carboni and Nicolas et al.
            %   (2017), Imaging of neuronal tissues by x-ray diffraction 
            %   and x-ray fluorescence microscopy: evaluation of contrast 
            %   and biomarkers for neurodegenerative diseases, Biomedical
            %   Optics Express 8 (10), pp. 4331-4347.
            %
            %       PARAMETER ARGUMENTS: 'symmetryArgs',struct()
            %           struct() can take the following arguments:
            %           - Description missing
            %           -
            %           -
            %           -         
            %
            %   heal: Heals a scattering pattern by filling values were
            %   data has been lost or not recorded by reflecting measured
            %   data from the point-symmetric location.
            %   Note, that heal is executed before all other methods, such
            %   that the healed data can be used for the subsequent analysis.
            %
            %       PARAMETER ARGUMENTS: 'healMask': 2d logical masks that
            %       identifies bad or corrupted pixels that should be
            %       recovered.
            % 
            %
            %   fluo: Integrates response within a pre-defined channel 
            %   range for each given counter. Note, that fluoCounters have
            %   to be defined beforehand. Either during initialization or
            %   afterwards using e.g. > obj.fluoCounters = {{'caka'},{200
            %   500}};
            %
            %   average: Averages a sequence of scattering data.
            %
            %   sum: Takes the sum of a sequence of scattering data.
            %
            %   data:  Simply outputs the data in a cell array. 
            %
            %   Scanning mode:
            %    Several raster scanning modes are possible:
            %     1 horizontal <horz> (default)
            %        *===================>
            %        ====================>
            %        ========>...            
            %     2 horizontal, alternating <horzalt>
            %        *===================>
            %                            v
            %        <===================<
            %        v
            %     3 vertical <vert>
            %        * # # #
            %        # # # #
            %        # # # #
            %        v v v v
            %     4 vertical,alternating <vertalt>
            %        * A>v A>
            %        # # # #
            %        # # # #
            %        v>A v>A
            %
            %   Output arguments:
            %       Result:: 
            %           struct containing the results for each method.
            %           Syntax: struct_name.method_name.method_results
            
            % make checks on input
            checkroi = @(x) isvector(x) && length(x) == 2;
            
            % read in data 
            pin = inputParser;
            % scan parameters
            addOptional(pin,'roi',[]);
            addOptional(pin,'yCrop',[1 obj.scan.SNy], checkroi);
            addOptional(pin,'zCrop',[1 obj.scan.SNz], checkroi);
            addOptional(pin,'ySkip',1,@isnumeric);
            addOptional(pin,'zSkip',1,@isnumeric);
            addOptional(pin,'scanmode','horz',@ischar);
            % general parameters
            addOptional(pin,'method','none',@ischar);
            addOptional(pin,'emptySub','off',@ischar);
            addOptional(pin,'empty',[]);
            addOptional(pin,'correction','off',@ischar);
            addOptional(pin,'parallel','off',@ischar);
            % stxm, crystal, b1d and symmetry
            addOptional(pin,'stxmArgs',struct(),@isstruct);
            addOptional(pin,'crystalArgs',struct(),@isstruct);
            addOptional(pin,'symmetryArgs',struct(),@isstruct);
            addOptional(pin,'b1dArgs',struct(),@isstruct);
            addOptional(pin,'pcaArgs',struct(),@isstruct);
            addOptional(pin,'pyfaiArgs',struct(),@isstruct);
            % fluo
            addOptional(pin,'fluoCounters',{'total',60,2300},@iscell);
            % healing
            addOptional(pin,'healMask',obj.mask);

            parse(pin,varargin{:});

            % shorthand notation
            zCrop = pin.Results.zCrop;
            yCrop = pin.Results.yCrop;
            ySkip = pin.Results.ySkip;
            zSkip = pin.Results.zSkip;
            method = pin.Results.method;
            scanmode = pin.Results.scanmode;
            
            % Process empty and masks
            emptyYes = 0;
            if strcmp(pin.Results.emptySub,'on')
                emptyYes = 1;
                if ~isempty(pin.Results.empty)
                    empty = pin.Results.empty;
                else
                    empty = (obj.empty);
                end
                
                if min(size(empty) ~= [obj.Nz,obj.Ny])
                    warning('Empty size does not match detector size');
                    try
                        empty = obj.process(empty);
                    catch err
                        error(err);
                    end
                end
            end
            if min(size(obj.mask) ~= [obj.Nz,obj.Ny])
                warning('Detector mask size does not match detector size');
                try
                    mask = obj.process_mask(obj.mask);
                catch err
                    error(err);
                end
            else
                mask = logical(obj.mask);
            end
            
            if min(size(obj.corr) ~= [obj.Nz,obj.Ny]) 
                warning('Correction mask size does not match detector size');
                try
                    corr = obj.process_mask(obj.corr);
                catch err
                    error(err);
                end
            else
                corr = obj.corr;
            end
            
            if ~isempty(pin.Results.roi)
                disp('Ignoring yCrop, zCrop, ySkip and zSkip since a custom roi is chosen');
                roi = pin.Results.roi;
                yCrop = [1 obj.scan.SNy];
                zCrop = [1 obj.scan.SNz];
                ySkip = 1;
                zSkip = 1;
            else
                % an automatic roi will be generated
                roi = zeros(obj.scan.SNz, obj.scan.SNy);
                roi(zCrop(1):zSkip:zCrop(2),yCrop(1):ySkip:yCrop(2)) = 1;                
            end

            % sampling points
            zsampling = zCrop(1):zSkip:zCrop(2);
            ysampling = yCrop(1):ySkip:yCrop(2);
            
            % Generate file numbers from roi                
            % get list of image numbers, starting at (1,1) with index 1.
            % Note, that linear indexing proceeds column-wise
            switch scanmode
                case 'horz'
                    indexlist = find(transpose(roi));
                case 'horzalt'
                    roi(2:2:end,:) = fliplr(roi(2:2:end,:));
                    indexlist = find(transpose(roi)); % not correct yet
                case 'vert'
                    indexlist = find(roi);
                case 'vertalt'
                    roi(:,2:2:end) = flipud(roi(:,2:2:end));
                    indexlist = find(roi); % not correct yet
                otherwise
                    error('Choose between <horz>,<horzalt>,<vert>,<vertalt>');
            end
            
            % Initialize Switches
            dataYes = 0; 
            rotateYes = 0; 
            pcaYes = 0; 
            stxmYes = 0; 
            develYes = 0;
            crystalYes = 0; 
            fluoYes = 0; 
            maxProjYes = 0; 
            averageYes = 0; 
            sumYes = 0; 
            symmetryYes = 0;
            b1dYes = 0;
            pyfaiYes = 0;
            healYes = 0;
            dpcYes = 0;
            liveYes = 0;
            
            nr_rows = length(zsampling);
            nr_cols = length(ysampling);
            if contains(method,'live')
                liveYes = 1;
            end
            if contains(method,'dpc')
                dpcYes = 1;
                dpcx = zeros(nr_rows*nr_cols,1);
                dpcy = zeros(nr_rows*nr_cols,1);
            end
            if contains(method,'heal')
                healYes = 1;
                healMask = pin.Results.healMask;
                [~, healIl] = heal(obj, zeros(obj.Nz,obj.Ny), healMask);
            end
            if contains(method,'devel')
                develYes = 1;
                devel = zeros(nr_rows*nr_cols,1);
            end
            if contains(method,'rotate') 
                rotateYes = 1;
            end
            if contains(method,'pca')
                pcaYes = 1;
                angle = zeros(nr_rows*nr_cols,size(obj.sels,3));
                w = zeros(nr_rows*nr_cols,size(obj.sels,3));
                pcaArgs = update_defaults(...
                    struct('qy',repmat(obj.qyaxis ,obj.Nz,1),'qz',repmat(obj.qzaxis',1,obj.Ny),'sel',obj.sels,'mask',obj.mask),...
                    pin.Results.pcaArgs);
            end
            if contains(method,'stxm') 
                stxmYes = 1;
                stxmArgs = update_defaults(...
                    struct('mask',obj.mask,'qrMin',-1,'qrMax',-1,'sel',ones(obj.Nz,obj.Ny)),...
                    pin.Results.stxmArgs);
                if (stxmArgs.qrMin > 0 && stxmArgs.qrMax < 0) 
                    stxmArgs.sel = obj.radial_mask('grid',obj.qr,'r1',stxmArgs.qrMin,'r2',1000);
                elseif (stxmArgs.qrMin < 0 && stxmArgs.qrMax > 0) 
                    stxmArgs.sel = obj.radial_mask('grid',obj.qr,'r1',0,'r2',stxmArgs.qrMax);
                elseif (stxmArgs.qrMin > 0 && stxmArgs.qrMax > 0) 
                    stxmArgs.sel = obj.radial_mask('grid',obj.qr,'r1',stxmArgs.qrMin,'r2',stxmArgs.qrMax);
                else
                end
                df = zeros(nr_rows*nr_cols,1);
            end
            if contains(method,'maxProj') 
                maxProjYes = 1;
                maxproj = zeros(obj.Nz,obj.Ny);
            end
            if contains(method,'crystal')
                crystalYes = 1;
                crystal = zeros(nr_rows*nr_cols,1);
                crystalArgs = update_defaults(...
                    struct('bgr',zeros(obj.Nz,obj.Ny),'mask',obj.mask,'threshold',10,'qrMin',0,'qrMax',100,'sel',ones(obj.Nz,obj.Ny)),...
                    pin.Results.crystalArgs);
                if (crystalArgs.qrMin > 0 && crystalArgs.qrMax < 0)
                    crystalArgs.sel = obj.radial_mask('grid',obj.qr,'r1',crystalArgs.qrMin,'r2',1000);
                elseif (crystalArgs.qrMin < 0 && crystalArgs.qrMax > 0)
                    crystalArgs.sel = obj.radial_mask('grid',obj.qr,'r1',0,'r2',crystalArgs.qrMax);
                elseif (crystalArgs.qrMin > 0 && crystalArgs.qrMax > 0)
                    crystalArgs.sel = obj.radial_mask('grid',obj.qr,'r1',crystalArgs.qrMin,'r2',crystalArgs.qrMax);
                else
                end
            end
            if contains(method,'average') 
                averageYes = 1; 
                sumYes = 1; 
                sumResult = zeros(obj.Nz,obj.Ny);
            end
            if contains(method,'sum') 
                sumYes = 1;
                sumResult = zeros(obj.Nz,obj.Ny);
            end
            if contains(method,'data') 
                dataYes = 1;
                dataStorage = cell(nr_rows * nr_cols, 1);
            end
            if contains(method,'symmetry') 
                symmetryYes = 1;
                Ln = zeros(nr_rows,nr_cols,length(0:4:ceil(max(obj.R(:)))));               
                symmetryArgs = update_defaults(...
                    struct('bgr',zeros(obj.Nz,obj.Ny),'sigma',1.6,'noise_level',6),...
                    pin.Results.symmetryArgs);
                % (optional) delta_theta = 20;
                % (optional) peak_thresh = 3;
            end
            if contains(method,'b1d') 
                b1dYes = 1;
                b1d_res = repmat(struct(), nr_rows * nr_cols, 1);
                b1dArgs = update_defaults(...
                    struct('grid',obj.qr,'bins',100,'mask',obj.mask,'sel',ones(obj.Nz,obj.Ny)),...
                    pin.Results.b1dArgs);
            end
            if contains(method,'fluo') 
                fluoYes = 1;
                obj.data.p.fluoFramesPerFile = obj.scan.SNy;
                % initialize struct fields for results corresponding to
                fluoResult = struct;
                fluoResult.averageSpectrum = zeros(4096,1);
                for ii = 1:length(obj.fluoCounters)
                    tmp = obj.fluoCounters{ii}(1);
                    % register field
                    fluoResult.(tmp{1}) = zeros(nr_rows,nr_cols);
                end
            end
            if contains(method,'pyfai') 
                pyfaiYes = 1;
                % Do not forget to import modules using init_pyfai
                pyfai_res = repmat(struct(), nr_rows * nr_cols, 1);
                
                pfsettings = update_defaults(...
                    struct('bins',100,'mask',[],'cake','off',...
                           'angles',[],'to2D','off','roi2D',[],...
                           'cakeRoi',[0 360]),...
                    pin.Results.pyfaiArgs);

                % cake integration
                if strcmp(pfsettings.cake,'on')
                    phi1 = pfsettings.cakeRoi(1);
                    phi2 = pfsettings.cakeRoi(2);
                    cake = py.tuple({phi1,phi2});
                    pfsettings.phi1 = phi1;
                    pfsettings.phi2 = phi2;

                    kwa_model = pyargs(...
                        'correctSolidAngle',0,...
                        'method','numpy',...
                        'unit','r_mm',...
                        'error_model','poisson',...
                        'azimuth_range',cake);  

                    pfsettings.kwa_model = kwa_model; 
                end                                        
            end
            
            
            tic;       
            
            if ~strcmp(pin.Results.parallel,'on')
                % serial version
                fprintf(1,'%6.2f%%', 0.00);
                for jj = 1:length(indexlist)

                    % progress 
                    fprintf(1,'\b\b\b\b\b\b\b%6.2f%%', double(jj)/length(indexlist)*100);

                    index = indexlist(jj);
                    fn = obj.scan.fnrs(index);
                    
                    % read data if necessary/requested
                    if ~contains(method,'none')

                        dat = obj.process(obj.data.read(fn));
%                         dat = obj.read_wait(fn);
                        
                        if emptyYes
                            dat = dat - empty;
                        end
                        
%                         dat = (~obj.mask).*data_temp;
                        
                        if strcmp(pin.Results.correction,'on')
                            dat = dat.*corr;
                        end
                        
                        % always apply mask
                        dat(mask) = 0;
                    end
                    
                    if healYes
%                         dat = obj.make_data_symmetric((~obj.mask).*dat,healMask);
                        dat = dat(healIl);
                    end
                    if liveYes
                        figure(1);imagesc(dat);colorbar;axis image;caxis([0 20]);
                        drawnow;
                    end
                    if dataYes
                        dataStorage{jj} = dat;
%                         dataStorage{jj} = (~obj.mask).*dat;
                    end
                    if fluoYes
                        fluoData = obj.data.read_fluorescence(index);
                        for kk = 1:length(obj.fluoCounters)
                            currCounter = obj.fluoCounters{kk}(1);
                            lower = obj.fluoCounters{kk}(2);
                            upper = obj.fluoCounters{kk}(3);
                            fluoResult.(currCounter{1})(index) = sum(fluoData(lower{1}:upper{1},1));
                        end
                        fluoResult.sumSpectrum = fluoResult.sumSpectrum + fluoData;
                    end
                    if maxProjYes
                        maxproj = max(dat,maxproj);
                    end
                    if stxmYes
                        for ii = 1:size(stxmArgs.sel,3)
                            df(jj,ii) = stxm(dat,stxmArgs.mask,stxmArgs.sel(:,:,ii));
                            % also available as mex function
                        end
                    end
                    if dpcYes
                        [dpcx(index),dpcy(index)] = dpc((~obj.mask).*dat);
                    end
                    if develYes
                        devel(index) = sum(sum(((~obj.mask).*dat) > 50));
                    end
                    if crystalYes
                        crystal(index) = crystals(dat-crystalArgs.bgr,crystalArgs.mask,crystalArgs.sel,crystalArgs.threshold);
                    end
                    if pcaYes
                        for ii = 1:size(pcaArgs.sel,3)
                            res = pca(dat,pcaArgs.qy,pcaArgs.qz,pcaArgs.mask,pcaArgs.sel(:,:,ii));
                            w(jj,ii) = res.w;
                            angle(jj,ii) = res.angle;
                            % also available as mex function
                        end
                    end
                    if sumYes
                        sumResult = sumResult + dat;
                    end
                    if symmetryYes
                        Ln(index,:) = symmetry(dat, symmetryArgs.bgr, symmetryArgs.sigma, symmetryArgs.noise_level, obj);
                    end
                    if b1dYes
                        for ii = 1:size(b1dArgs.sel,3)
                            tmp = b1d(dat,b1dArgs.mask,b1dArgs.sel(:,:,ii),b1dArgs.grid,b1dArgs.bins);
                            b1d_res(index,ii).dat_1d = tmp.dat_1d;
                            b1d_res(index,ii).qr = tmp.qr;
                        end
                    end
                    if pyfaiYes
                        tmp = pyfai(dat,pin.Results.pyFAImask,obj.pby,obj.pbz,pfsettings);
                        pyfai_res(index).error = tmp.error;
                        pyfai_res(index).dat_1d = tmp.dat_1d;
                        pyfai_res(index).qr = obj.q_of_n(tmp.r);
                        if strcmp(pfsettings.to2D,'on')
                            pyfai_res(index).dat_2d = tmp.dat_2d;
                        end
                    end
                        
                end
                fprintf(1,'\b\b\b\b\b\b\b');
            else 
                % here parallel version
                
                f = obj.data.compile_static_function;
                g = obj.process_data_static;
                
                % slice some variables
                emptySub = pin.Results.emptySub;
                fnrs = obj.scan.fnrs;
                mask_sliced = (~obj.mask);
                parfor jj = 1:length(indexlist)

                    index = indexlist(jj);
                    fn = fnrs(index);
                 
                    % always read data
                    data_temp = f(fn);

                    if strcmp(emptySub,'on')
                        dat = mask_sliced.*(g(data_temp) - empty);
                    else                            
                        dat = mask_sliced.*g(data_temp);
                    end
                    
                    if stxmYes
                        df(jj) = sum(dat(:));
                    end
                    
                end                
            end
            
            toc;
           
            % save result
            if dataYes
                result.data = sort_scan(dataStorage,nr_rows,nr_cols,scanmode);
            end
            if maxProjYes
                result.maxProjection = maxproj;
            end
            if crystalYes
                result.crystal = sort_scan(crystal,nr_rows,nr_cols,scanmode);
                result.crystalBgr = crystalBgr;
            end
            if develYes
                result.devel = devel;
            end
            if dpcYes
                result.dpc.dpcx = sort_scan(dpcx,nr_rows,nr_cols,scanmode);
                result.dpc.dpcy = sort_scan(dpcy,nr_rows,nr_cols,scanmode);
            end
            if stxmYes
                result.stxm.df = sort_scan(df,nr_rows,nr_cols,scanmode);
            end
            if pcaYes
                result.pca.w = sort_scan(w,nr_rows,nr_cols,scanmode);
                result.pca.angle = sort_scan(angle,nr_rows,nr_cols,scanmode);
            end
            if fluoYes
                result.fluo = fluoResult;
                for kk = 1:length(obj.fluoCounters)
                    currCounter = obj.fluoCounters{kk}(1);
                    fluoResult.(currCounter{1}) = sort_scan(fluoResult.(currCounter{1}),nr_rows,nr_cols,scanmode);
                end
            end
            if averageYes
                result.avg = sumResult / length(indexlist);
            end
            if sumYes
                result.sum = sumResult;
            end
            if rotateYes
                result.rotate = angles;
            end
            if symmetryYes
                result.symmetry.Ln = Ln;
            end
            if b1dYes
               result.b1d = sort_scan(b1d_res,nr_rows,nr_cols,scanmode);
            end
            if pyfaiYes
               result.pyfai = sort_scan(pyfai_res,nr_rows,nr_cols,scanmode);
            end
        end        
        
        
        function [healed, il] = heal(obj, data, mask)
            % HEAL  uses make_mask_symmetric to find the correct mapping
            % from valid to invalid pixels. Takes in an n x m data matrix
            % and an n x m logical mask (1: corrupt pixel, 0: otherwise)
            % and mirrors data that is lost due to e.g. modular gaps from 
            % the point-symmetric position in the data.
            %   
            %   HEALED = HEAL(DATA,MASK)
            %
            %   The following parameters are accepted:
            %
            %       data:: [] 
            %          An n x m data matrix.
            %                
            %       mask:: [] 
            %          An n x m logical matrix indicating bad pixels.
            %          1: invalid pixel, 0: valid pixel.
            %
            %   Output arguments:
            %   
            %       healed::
            %           healed scattering pattern.
            %
            %       il:: 
            %           index list (mapping from valid to invalid pixels)
            
            % indexlist
            il = transpose(1:obj.Ny*obj.Nz);
            il = reshape(il,obj.Nz,obj.Ny);
            
            [~,~,~,good_zero,bad_zero] = obj.make_mask_symmetric(mask);
            bad_zero = logical(bad_zero);
            good_zero = logical(good_zero);
            
            il(bad_zero) = flipud(il(good_zero));
            healed = data(il);
        end
        
        
        function data = read_wait(obj,fn)
            % READ_WAIT  Useful for online data analysis when not all data
            % is immediately ready for processing. Similar to read() but
            % uses (nanodiffraction).wait_it and
            % (nanodiffraction).wait_maxit and (nanodiffraction).wait to
            % decide wether or not to wait (wait) and for how long 
            % (wait_maxit). wait_it is the current waiting index.
            %   
            %   DATA = READ_WAIT(FN)
            %
            %   Input and output arguments as in read().
            
            
            % should I wait?
            if obj.wait && (obj.wait_it < obj.wait_maxit)
                try
                    data = obj.process(obj.data.read(fn)); 
                    obj.wait_it = 0;
                catch err
                    obj.wait_it = obj.wait_it + 1;
                    pause(1);
                    data = obj.read(fn);
                end 
            else
                data = obj.process(obj.data.read(fn)); 
            end
        end
        
        
        function [composite, parameters] = calculate_composite(obj,varargin)
            % CALCULATE_COMPOSITE  calculates a composite image from a
            % series of diffraction patterns.
            %   
            %   [COMPOSITE, PARAMETERS] = CALCULATE_COMPOSITE(OPTS) 
            %
            %   The following options are supported:
            %
            %     OPTS:: [] (optional)
            %       Struct that contains the parameters that were used for 
            %       composite calculation.
            %       p can be obtained as a second argument from
            %       e.calculate_composite() or can be defined manually.
            %       DEFAULT ARGUMENTS: 
            %           'yCrop': [1 obj.scan.SNy]
            %           'zCrop': [1 obj.scan.SNz]
            %           'ySkip': 1
            %           'zSkip': 1
            %
            %   Output arguments:
            %       COMPOSITE::
            %           Processed composite matrix.
            %
            %       PARAMETERS::
            %           Parameters used to combine the input data. Is used
            %           in combination with (display_class).composite. See
            %           help display.composite for more information
                        
            defaults = struct('yCrop',[1 obj.scan.SNy],'zCrop',[1 obj.scan.SNz],'ySkip',1,'zSkip',1,'correction','off','emptySub','off','empty',[],'heal','off','healmask',[]);
            fields = fieldnames(defaults);
            if nargin > 1
                opts = varargin{1};
                for f = 1:numel(fields)
                    if ~isfield(opts,fields{f})
                        opts.(fields{f}) = defaults.(fields{f});
                    end
                end
            else
                opts = defaults;
            end
            
            % size of the detector image
            Ny = obj.Ny;
            Nz = obj.Nz;

            % get data stack
%             if strcmp(opts.heal,'off')
%                 result = obj.analyze_scan('method','data','yCrop',opts.yCrop,'zCrop',opts.zCrop,'ySkip',opts.ySkip,'zSkip',opts.zSkip,'emptySub',opts.emptySub,'correction',opts.correction);
%             else
                result = obj.analyze_scan('method','data+heal','yCrop',opts.yCrop,'zCrop',opts.zCrop,'ySkip',opts.ySkip,'zSkip',opts.zSkip,'emptySub',opts.emptySub,'empty',opts.empty,'correction',opts.correction,'healmask',opts.healmask);
%             end
            
            % yCrop and zCrop, ySkip and zSkip are the important parameters
            imsY = numel(opts.yCrop(1):opts.ySkip:opts.yCrop(2));
            imsZ = numel(opts.zCrop(1):opts.zSkip:opts.zCrop(2));
            
            % final matrix
            composite = zeros(imsZ*Nz,imsY*Ny);
            for r = 1:imsZ
                for c = 1:imsY
                    composite((1:Nz) + (r-1)*Nz,(1:Ny) + (c-1)*Ny) = result.data{r,c};
                end
            end
            
            parameters = opts;
            parameters.Ny = Ny;
            parameters.Nz = Nz;
            parameters.imsY = imsY;
            parameters.imsZ = imsZ;
            parameters.binning = obj.binning;
            parameters.detRoiY = obj.detRoiY;
            parameters.detRoiZ = obj.detRoiZ;
            
        end
            
        function [bgr] = median_background(obj)
            % MEDIAN_BACKGROUND  Calculates the median of the first lines
            % in a scan. In some cases, this can be used as a simple
            % estimate for a background for scattering patterns.
            %   
            %   [BACKGROUND] = MEDIAN_BACKGROUND() 
            %
            %   Output arguments:
            %       BACKGROUND::
            %           Calculated median background.
            
            
            % line median
            fprintf(1,'Averaging %d lines for median calculation...', ceil(50/double(obj.scan.SNy)));
            rawdata = obj.analyze_scan('zCrop',[1 ceil(50/double(obj.scan.SNy))],'method','data');

            % reshape
            [m,n] = size(cell2mat(rawdata.data(1)));
            stack = zeros(m,n,numel(rawdata.data));
            for ii = 1:numel(rawdata.data)
                stack(:,:,ii) = cell2mat(rawdata.data(ii));
            end

            % median-filtered background
            bgr = median(stack,3);
        end
        
        
        function robustFit = robust_fitting(obj,dat,varargin)
            % ROBUST_FITTING  in principle, this function performs an 
            % ordinary weighted linear least square analysis. Since it is
            % linear, it can be written in a matrix equation form. 
            % For more information, see:
            % https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/06LeastSquares/general/
            % http://de.mathworks.com/help/curvefit/least-squares-fitting.html
            % http://www.dsplog.com/2012/02/05/weighted-least-squares-and-locally-weighted-linear-regression/
            % 
            %   [ROBUST_FIT] = ROBUST_FITTING(DATA,OPTS)
            %
            %   The following options are supported:
            %
            %     DATA:: []
            %       A struct containing the fields 'dat_1d' and 'qr' is
            %       required.
            %
            %     OPTS:: [] (optional)
            %       structure that can contain the following fields:
            %           - 'ncoeff': number of coefficients
            %           - 'win': Window for data selection, a 2-element vector.
            %
            %   Output arguments:
            %    
            %       ROBUST_FIT::
            %           struct containing
            %           - f: fit amplitude
            %           - x: fit abscissa
            %           - c: fitted coefficients
            
            % parse variable arguments
            defaults = struct('ncoeff',9,'win',[1 numel(dat.qr)],'coeff_start',0);
            fields = fieldnames(defaults);
            if nargin > 2
                opts = varargin{1};
                for f = 1:numel(fields)
                    if ~isfield(opts,fields{f})
                        opts.(fields{f}) = defaults.(fields{f});
                    end
                end
            else
                opts = defaults;
            end
            
            n_coeff = opts.ncoeff;
            win = opts.win;

            % data
            x = dat.qr(win(1):win(2));
            y = dat.dat_1d(win(1):win(2));
            
            % fit
            [f,c] = lin_fit(x,y, n_coeff, @(y) 1./y.^2, opts.coeff_start);

            % save
            robustFit.f = f;
            robustFit.x = x;
            robustFit.c = c;
        end
        
        
        function [res] = cluster(~,c,max_clusters,varargin)
            % CLUSTER  This function starts with k random starting points 
            % in the nx x ny x c_i - space 
            % (e.g. a scan with 51x51 scan points and 9 fitting
            % coefficients gives a 51x51x9 volume). It then searches for
            % k means within this volume by minimizing the euklidian norm.
            % 
            %   [CLUSTERS] = CLUSTER(C,MAX_CLUSTERS,LOGSCALE)
            %
            %   The following options are supported:
            %
            %     C:: []
            %       An n x m x k - data matrix where k = 1:K and K is the 
            %       number of clusters. n and m are row and col position in
            %       the scan.
            %
            %     MAX_CLUSTERS:: [] 
            %       Maximum numbers of clusters (i.e. means).
            %
            %     LOGSCALE:: [1] (optional)
            %       If set to 1, k-means will be displayed in log-scale.
            %
            %   Output arguments:
            %    
            %       See help cluster_analysis.
            
            if nargin == 4
                logscale = varargin{1};
            else
                logscale = 1;
            end
            
            nx = size(c,2);
            ny = size(c,1);
            nz = size(c,3);
            
            % bring c into form of (coeff,datapoint)
            c_reshaped = zeros(nz,nx*ny);
            for ii = 1:size(c,3)
                c_reshaped(ii,:) = squeeze(reshape(c(:,:,ii),[1,nx*ny]));
            end
            
            res = cluster_analysis(c_reshaped,max_clusters,logscale);
            res.clusters = reshape(res.cluster,ny,nx);
            
        end
        
        function mask = process_mask(obj,mask)
            % PROCESS_MASK  Used for processing logical masks. See help
            % nanodiffraction.mask for usage.
            
            mask = obj.process(mask);
            mask = min(mask,1);
        end
            
        function selections = add_selection(obj, sel)
            % ADD_MASK  The function accepts a task and a mask. The mask
            % will be added to the stack of masks corresponding to the
            % given task.
            % 
            %   [SELECTIONS] = ADD_SELECTION(SEL)
            %
            %   The following options are supported:
            %
            %     SEL:: []
            %       A two-dimensional array. The input should be a logical
            %       array. 
            %
            %   Output arguments:
            %
            %       SELECTIONS::
            %           Current set of selections.
            
            % add selection if dimensions fit to previously added selection
            if numel(sel) == numel(obj.sels(:,:,1))
                obj.sels = cat(3,obj.sels,sel);
            else
                error('Selection matrix does not have the same size as the selections stored in obj.sels');
            end

            selections = obj.sels;
            
        end
        
                    
        function mask_out = set_mask(obj, varargin)
            % SET_MASK  The function accepts as many masks as needed as
            % arguments, as long as their dimensions are identical. 
            % Note, that a warning will be shown, if the size of the mask 
            % does not correspond to the current detector dimensions. 
            % 
            %   [MASK_OUT] = SET_MASK(TASK, MASK, MASK, MASK, ...)
            %
            %   The following options are supported:
            %
            %     MASK:: [] (required)
            %       A two-dimensional array.The input should be a logical
            %       array. 
            %
            %     MASK:: [] (optional)
            %       Additional masks can be given. All masks will then be 
            %       logically combined (OR combination [|]).
            %
            %     MASK:: [] (optional)
            %       The amount of masks is not limited, all masks will be
            %       combined to one.
            %
            %   Output arguments:
            %
            %       MASK_OUT::
            %           Combined mask.
            %
            %   Note, that the detector mask is stored in:
            %       (nanodiffraction_class).mask
            
            % combination of masks, only applicable on logical masks
            mask = zeros(size(varargin{1}));
            for ii = 1:(nargin-1)
                mask = mask + varargin{ii};
            end
            mask = min(mask,1); % clip to 1
            
                
            if max(size(mask) ~= [obj.Nz,obj.Ny])
                warning('Size of detector mask and dimensions of the detector do not match. Mask will be processed to match the current roi and binning settings.');
                mask = obj.process_mask(mask);
            end
            
            obj.mask = mask;
            
            if max(size(mask) == [obj.Nz_orig,obj.Ny_orig])
                obj.mask_orig = obj.mask; 
            end
            
            mask_out = mask;
        end
        
        function sel_out = set_selection(obj, varargin)
            % SET_SELECTION  The function accepts as many masks as needed as
            % arguments, as long as their dimensions are identical. 
            % Note, that a warning will be shown, if the size of the mask 
            % does not correspond to the current detector dimensions. 
            % 
            %   [SEL_OUT] = SET_SELECTION(SEL, SEL, SEL, ...)
            %
            %   The following options are supported:
            %
            %     SEL:: [] (required)
            %       A two-dimensional array. The input should be a logical 
            %       array. 
            %
            %     SEL:: [] (optional)
            %       Additional selections can be given. All selections will 
            %       then be logically combined (OR combination [|]).
            %
            %     SEL:: [] (optional)
            %       The amount of selections is not limited, all selections 
            %       will be combined to one.
            %
            %   Output arguments:
            %
            %       SEL_OUT::
            %           Current selection.
            %
            %   Note, that selections are stored in:
            %       (nanodiffraction_class).sels
            %
            %   Note, that more selections can be used, for example, when
            %   performing a pca analysis or stxm. Additional selections
            %   can be added by calling 
            %       (nanodiffraction).add_selection(SEL)
            %
            %   For more help, please read:
            %       help nanodiffraction.add_selection
            
            % combined selections, only applicable on logical masks
            sel = zeros(size(varargin{1}));
            for ii = 1:(nargin-1)
                sel = sel + varargin{ii};
            end
            sel = min(sel,1); % clip to 1
            
            if max(size(sel) ~= [obj.Nz,obj.Ny])
                warning('Size of detector mask and dimensions of the detector do not match. Mask will be processed to match the current roi and binning settings.');
                sel = obj.process_mask(sel);
            end
            
            obj.sels = sel;     
            
            if max(size(sel) == [obj.Nz_orig,obj.Ny_orig])
                obj.sels_orig = obj.sels;     
            end
            
            sel_out = sel;
        end
        
        function corr_out = set_corr(obj, corr)
            % SET_CORR  Defines a correction matrix to be applied in a
            % subsequent analysis.
            % Note, that a warning will be shown, if the size of the
            % correction matrix does not correspond to the current detector 
            % dimensions. 
            % 
            %   [CORR_OUT] = SET_CORR(CORR)
            %
            %   The following options are supported:
            %
            %     CORR:: [] 
            %       A two-dimensional array. Input can a double array.
            %
            %   Output arguments:
            %       CORR_OUT::
            %           Active correction matrix
            %
            %   Note, that the correction matrix is stored in: 
            %   (nanodiffraction_class).corr
            %   after set_corr was used.
            
 
            if max(size(corr) ~= [obj.Nz,obj.Ny])
                warning('Size of correction array and dimensions of the detector do not match. Correction array will be processed to match the current roi and binning settings.');
                corr = obj.process_data(corr);
            end
            
            corr(isinf(corr)) = 0; % division by 0 results in Infs
            corr(isnan(corr)) = 0; % interpolated, might contain NaNs
            obj.corr = corr;
                            
            if max(size(corr) == [obj.Nz_orig,obj.Ny_orig])
                obj.corr_orig = obj.corr;
            end
            
            corr_out = corr;
        end
        
       
        function obj = set_roi_and_binning(obj,varargin)
            % SET_DETECTOR_MASK  All following functions require a standard 
            % pixel mask that mask hot pixels and intermodular gaps of 
            % detectors. In addition, a calibration could be added that 
            % corrects for solid angle deviations at larger scattering 
            % angles.
            %
            %   The following arguments are accepted:
            %
            %      detCalibration:: [1] 
            %        mask taking into account corrections e.g. for solid 
            %        angles
            %
            %      detectorRoi:: [off]
            %        Can be toggled between on and off
            %
            %      detRoiY:: [start end] in pixel units
            %        [start end] in pixel units
            %
            %      detRoiZ:: [start end] in pixel units
            %        [start end] in pixel units  
            %
            %      binning:: [off]
            %        Can be toggled between on and off
            %
            %      binx:: [1]
            %        Binning ratio along y
            %
            %      biny:: [1]
            %        Binning ratio along z
        
            % parse input
            pin = inputParser;
            addOptional(pin,'detectorRoi','off');
            addOptional(pin,'detRoiY',[1 obj.Ny]);            
            addOptional(pin,'detRoiZ',[1 obj.Nz]);
            addOptional(pin,'binning','off');
            addOptional(pin,'biny',1);
            addOptional(pin,'binz',1);
            parse(pin,varargin{:});
            
            % reset the values to their original value when the class was
            % initialized
            obj.Ny = obj.Ny_orig;
            obj.Nz = obj.Nz_orig;
            obj.pby = obj.pby_orig;
            obj.pbz = obj.pbz_orig;
            obj.pixelsize = obj.pixelsize_orig;
            
            % if roi or binning was already activated, masks are already
            % processed, therefore, if set_roi_and_binning is called a
            % second time, a warning should be issued and the original
            % matrices should be used.
            if strcmp(obj.detectorRoi,'on') || strcmp(obj.binning,'on')
                % either roi or binning is toggeled off
                if strcmp(pin.Results.detectorRoi,'off') || strcmp(pin.Results.binning,'off')
                    warning('Masks are being reset to their initial value');
                    obj.mask = obj.mask_orig;
                    obj.sels = obj.sels_orig;
                    obj.corr = obj.corr_orig;
                else
                    warning('Do you attempt to set another roi or binning ratio?');
                end
            end
            
            % save 
            obj.detectorRoi = pin.Results.detectorRoi;
            obj.detRoiY = pin.Results.detRoiY;
            obj.detRoiZ = pin.Results.detRoiZ;
            obj.binning = pin.Results.binning;
            obj.binx = pin.Results.biny;
            obj.biny = pin.Results.binz;

            if strcmp(obj.detectorRoi,'on')
                % make the according changes to... (detector Roi)
                obj.Ny = (obj.detRoiY(2) - obj.detRoiY(1)) + 1;
                obj.Nz = (obj.detRoiZ(2) - obj.detRoiZ(1)) + 1;
                obj.pby = obj.pby - obj.detRoiY(1) + 1;
                obj.pbz = obj.pbz - obj.detRoiZ(1) + 1;
            end
            
            if strcmp(obj.binning,'on')
                obj.Ny = (obj.Ny - mod(obj.Ny,obj.binx))/obj.binx;
                obj.Nz = (obj.Nz - mod(obj.Nz,obj.biny))/obj.biny;
                obj.pby = obj.pby/obj.binx;
                obj.pbz = obj.pbz/obj.biny;
                obj.pixelsize = obj.pixelsize.*obj.binx;
            else
                obj.binx = 1;
                obj.biny = 1;
            end
            
            % recalculate q_r and so on and so forth
            obj.refresh();
            
        end        
        
        
        function obj = set_scan_info(obj,varargin)
            % SET_SCAN_INFO  sets the scan information for all subsequent 
            % analyses.
            %
            %   The following arguments can be set:
            %   
            %       SNy, SNz:: [1]
            %          number of scan points along fast and slow axis
            %   
            %       stepy, stepz:: [1e-6]
            %          Stepsize along both axes in meters
            %
            %       firstFile:: [1]
            %          In some special cases, the first file of a scan is
            %          not 1. If this is the case, please change the
            %          "firstFile" argument accordingly (however, this 
            %          should not normally happen)
            %
            % For p10, one could use spec_get_scan_info, therefore choose
            % p10 for the beamline argument. If one chooses to do so,
            % please also provide a "scanNo","newfile" and "specpath"
            % argument, such that "spec_get_scan_info" can be used in turn.
            % However, this is usually too complicated and the beamline
            % argument can be disregarded completely. 
            
            % parse input
            pin = inputParser;
            addOptional(pin,'scanNo',1);
            addOptional(pin,'newfile','');
            addOptional(pin,'specpath','');
            addOptional(pin,'SNy',1);
            addOptional(pin,'SNz',1);
            addOptional(pin,'stepy',1E-6);
            addOptional(pin,'stepz',1E-6);            
            addOptional(pin,'roi',0);
            addOptional(pin,'roiMask',0);
            addOptional(pin,'beamline','none');
            addOptional(pin,'firstFile',1);
            parse(pin,varargin{:});
                
            % scanning axis
            y = (0:(pin.Results.SNy)-1) * pin.Results.stepy;
            z = (0:(pin.Results.SNz)-1) * pin.Results.stepz;

            % file numbers for each scanpoint
            fnrs = (1:1:(pin.Results.SNy*pin.Results.SNz)) + pin.Results.firstFile - 1;

            % save to object
            obj.scan.SNy = pin.Results.SNy;
            obj.scan.SNz = pin.Results.SNz;
            obj.scan.fnrs = fnrs;
            obj.scan.y = y;
            obj.scan.z = z;                          
            if (isfield(pin.Results,'roi') || isfield(pin.Results,'roiMask'))
                if isfield(pin.Results,'roiMask')
                   obj.scan.roiMask = pin.Results.roiMask;
                else
                   obj.scan.roi = pin.Results.roi;
                   obj.scan.roiMask = zeros(size(obj.scan.SNz,obj.scan.SNy));
                   obj.scan.roiMask(roi(3):roi(4),roi(1):roi(2)) = 1;
                end
            end

        end
        
        
        function empty = calculate_empty(obj,varargin)
            % CALCULATE_EMPTY  calculates an empty image and stores the
            % resulting empty in <experiment>.scan.empty
            %
            %   The following arguments are accepted:
            %
            %     roi:: []
            %        Logical mask that defines a region-of-interest
            %
            %      yCrop:: [1 NR_OF_SCAN_POINTS_ALONG_FAST_AXIS] 
            %        [start end] in pixel units
            %
            %      zCrop:: [1 NR_OF_SCAN_POINTS_ALONG_SLOW_AXIS] 
            %        [start end] in pixel units  
            %
            %      ySkip:: [1]
            %        only every n-th (n integer) scan point along fast axis
            %        shall be evaluated
            %
            %      zSkip:: [1]
            %        only every n-th (n integer) scan point along slow axis
            %        shall be evaluated
        
            % make checks on input
            checkroi = @(x) isvector(x) && length(x) == 2;
            
            % parse input
            pin = inputParser;
            addOptional(pin,'roi',[]);
            addOptional(pin,'yCrop',[1 obj.scan.SNy], checkroi);
            addOptional(pin,'zCrop',[1 obj.scan.SNz], checkroi);
            addOptional(pin,'ySkip',1,@isnumeric);
            addOptional(pin,'zSkip',1,@isnumeric);
            parse(pin,varargin{:});
            
            if isempty(pin.Results.roi)
                res = obj.analyze_scan('yCrop',pin.Results.yCrop,'zCrop',pin.Results.zCrop,'ySkip',pin.Results.ySkip,'zSkip',pin.Results.zSkip,'method','average');
            else
                res = obj.analyze_scan('roi',pin.Results.roi,'method','average');
            end
            
            % save the result
            obj.scan.empty = res.avg;
            empty = res.avg;
        end
        
        function set_empty(obj,empty)
            obj.scan.empty = empty;
        end
        
       
        function dat = process(obj,dat)
            % PROCESS_DATA  Post-processing of data if detector roi or
            % binning was toggeled to on.
            %
            %   DATA_OUT = PROCESS_DATA(DATA_IN)
            %
            % See set_detector_mask for more information
            
            
            if obj.detRoiZ(1) > size(dat,1) || obj.detRoiZ(2) > size(dat,1) || ...
                    obj.detRoiY(1) > size(dat,2) || obj.detRoiY(2) > size(dat,2)
                fprintf('Size of matrix: %d (row) x %d (col)\n', size(dat,1), size(dat,2));
                fprintf('Roi (y): [%d %d], Roi (z): [%d %d]\n', obj.detRoiY(1),obj.detRoiY(2),obj.detRoiZ(1),obj.detRoiZ(2));
                warning('Sizes do not match');
                return
            end
            
            % apply detectorRoi only if really desired
            if strcmp(obj.detectorRoi,'on')
                dat = dat(obj.detRoiZ(1):obj.detRoiZ(2),obj.detRoiY(1):obj.detRoiY(2));
            end
            
            % binning
            if strcmp(obj.binning,'on')        
                
                [m,n] = size(dat);
                nnew = n - mod(n,obj.binx);
                mnew = m - mod(m,obj.biny);

                % now it is divisible by binx, biny
                dat = dat(1:mnew,1:nnew);
                
                % binning along z (correct)
                dat=sum( reshape(dat,obj.biny,[]),1); 
                dat=reshape(dat,mnew/obj.biny,[]); 

                % binning along y
                dat=sum(reshape(dat',obj.binx,[]),1); 
                dat=reshape(dat,nnew/obj.binx,[])';
            end
        end
        
        
        function g = process_data_static(obj)
            % PROCESS_DATA_STATIC  slices all variables for use in parfor
            % loop (under development)
            %
            %   FUNCTIONHANDLE = PROCESS_DATA_STATIC()
            
            biny = obj.biny;
            binx = obj.binx;
            detRoiZ = obj.detRoiZ;
            detRoiY = obj.detRoiY;
            if strcmp(obj.detectorRoi,'on') && ~strcmp(obj.binning,'on') 
                g = @(dat) dat(detRoiZ(1):detRoiZ(2),detRoiY(1):detRoiY(2));
                return;
            end
            if strcmp(obj.binning,'on') && ~strcmp(obj.detectorRoi,'on')
                [m,n] = deal(obj.Nz_orig,obj.Ny_orig);
                nnew = n - mod(n,obj.binx);
                mnew = m - mod(m,obj.biny);
                g = @(dat) transpose(reshape(sum(reshape(transpose(reshape(sum( reshape(dat(1:mnew,1:nnew),biny,[]),1),mnew/biny,[])),binx,[]),1),nnew/binx,[]));
                return;
            end
            if strcmp(obj.binning,'on') && strcmp(obj.detectorRoi,'on')
%                 [m,n] = deal(obj.Nz,obj.Ny);
                [m,n] = deal(detRoiZ(2)-detRoiZ(1)+1,detRoiY(2)-detRoiY(1)+1);
                
                nn = n - mod(n,obj.binx);
                mm = m - mod(m,obj.biny);

                % now it is divisible by binx, biny
                g = @(dat) reshape(sum(sum(...
                    reshape(...
                        dat(detRoiZ(1):(detRoiZ(2) - mod(n,obj.binx)),detRoiY(1):(detRoiY(2) - mod(n,obj.binx))),...
                        [biny mm/biny binx nn/binx]),...
                    1),3),[mm/biny nn/binx]);

                return;
            else
                g = @(dat) dat;
                return;
            end
        end
        
        
        function mask = radial_mask(obj,varargin)
            % RADIAL_MASK  defines a radial mask (e.g. necessary for pca)
            %
            %   MASK = RADIAL_MASK(VARARGIN) creates a circular or ring-shaped
            %   mask and stores it in <experiment>.p.sels
            %
            %   The following parameters are accepted:
            %
            %       r1:: [65] 
            %          Inner radius
            %
            %       r2:: [232] 
            %          Outer radius (Inverse radius)
            %
            % If both parameters are given, the function will generate
            % logical 1, if r is within r1 and r2
            % If only r1 is given, logical mask will be 1 for r < r1
            % If only r2 is given, logical mask will be 1 for r > r2
        
            % parse input
            pin = inputParser;
            addOptional(pin,'grid',obj.R);
            addOptional(pin,'r1',65);
            addOptional(pin,'r2',232);
            parse(pin,varargin{:});
            
            % prepare radial mesh
            R = pin.Results.grid;
            mask = ones(size(R));
            
            % check which Defaults are used
            r1_provided = sum(cell2mat(strfind(fieldnames(pin.Results),'r1')));
            r2_provided = sum(cell2mat(strfind(fieldnames(pin.Results),'r2')));
            condition = (r1_provided + r2_provided);
            switch condition
                case 0
                    % no input arguments, throw warning that defaults will 
                    % be used
                    warning('No arguments provided, method will not be used')

                case 1    
                    if r2_provided
                        R_max = pin.Results.r2;
                        mask = ~(R - R_max <= 0);
                    end 
                    if r1_provided
                        R_max = pin.Results.r1;
                        mask = (R - R_max <= 0);
                    end
                    
                case 2
                    % no defaults given
                    if ~(pin.Results.r1 < pin.Results.r2)
                        error('r1 has to be smaller than r2')
                    end
                    R_max = pin.Results.r1;
                    mask1 = R - R_max <= 0;
                    R_max = pin.Results.r2;
                    mask2 = R - R_max <= 0;

                    % total mask
                    mask = ~mask1 .* mask2;
            end
        end
        
        
        function mask = azimuthal_mask(obj,varargin)
            % AZIMUTHAL_MASK  defines an azimuthal mask (e.g. necessary for pca)
            %
            %   MASK = AZIMUTHAL_MASK(VARARGIN) creates a circular or 
            %   ring-shaped mask. Note, that the angles in degrees follow 
            %   the mathematical convention (counter-clockwise and 0° at
            %   3 o'clock).
            %
            %   The following parameters are accepted:
            %
            %       phi1:: [65] 
            %          Start angle in degrees.
            %
            %       phi2:: [232] 
            %          End angle in degrees
            %
            % If both parameters are given, the function will generate
            % logical 1, if r is within r1 and r2
            % If only r1 is given, logical mask will be 1 for r < r1
            % If only r2 is given, logical mask will be 1 for r > r2
        
            % parse input
            pin = inputParser;
            addOptional(pin,'phi1',65);
            addOptional(pin,'phi2',232);
            parse(pin,varargin{:});
            
            % shortform
            phi1 = mod(pin.Results.phi1,360);
            phi2 = mod(pin.Results.phi2,360);
            
            % coordinate transform
            Y = repmat(obj.yaxis,obj.Nz,1);
            Z = repmat(obj.zaxis',1,obj.Ny);
            PHI = atan2d(Z,Y);
            PHI(PHI<0) = PHI(PHI<0) + 360;
            mask = ones(size(PHI));
            
            % check which Defaults are used
            phi1_provided = sum(cell2mat(strfind(fieldnames(pin.Results),'phi1')));
            phi2_provided = sum(cell2mat(strfind(fieldnames(pin.Results),'phi2')));
            condition = (phi1_provided + phi2_provided);
            if condition == 2
                    % no defaults given
                    if ~(phi1 < phi2)
                        mask1 = PHI > phi1;
                        mask2 = PHI < phi2;

                        % total mask
                        mask = min(mask1 + mask2,1);
                    else
                        mask1 = PHI < phi2;
                        mask2 = PHI > phi1;

                        % total mask
                        mask = mask1 .* mask2;
                    end
            else
                    % no input arguments, throw warning that defaults will 
                    % be used
                    warning('phi1 and phi2 have to be provided, method will not be used')
            end
        end
        
                
        function [cm] = estimate_primary_beam(obj,data)
            % ESTIMATE_PRIMARY_BEAM  function that estimates primary beam 
            % position based on center of mass (under development)
        
            % estimate primary beam position
            [obj.pbz,obj.pby] = CM_nano(data); % output format [col, row]

            % display result if needed
            prompt = 'Do you wish to see the result? (y/n)';
            answer = input(prompt,'s');
            if strcmp(answer,'y')
                imagesc(log10(abs(data)));
            end
        end
        
        
        function [orientation] = orientation(~,pca)
            % ORIENTATION  shorthand to get the orientation angle instead 
            % of the x- and y-component of the angle
            %
            %   PHI = ORIENTATION(PCA) accepts a struct that is formatted
            %   as the output of analyze scan, as 
            %       pca - output1
            %           - output2...
            
            orientation = atan(pca.ycomp./pca.xcomp);
        end
        
        
        function [circMean, circVar] = circular_dist(~,phi,I)
            % CIRCULAR_DIST  Calculates the circular mean and variance of a 
            % distribution based on 
            % https://en.wikipedia.org/wiki/Directional_statistics
            %
            % Accepted arguments:
            % phi, I, where I is usually given in counts and phi usually 
            % within the range [-180, 180] (because of central symmetry)
            
            angles = [];
            for jj = 1:numel(phi)
                % repeat the angle N times for each angle
                angles = [angles repmat(phi(jj),1,round(I(jj)))];
            end
            
            % number of data points
            N = numel(angles);
            
            % transform to complex number
            z = zeros(N,1);
            for jj = 1:N
                z(jj) = 1.0*exp(1i*angles(jj)*pi/180);
            end
            
            % circular mean
            circMean = angle(sum(z)/N)*180/pi;
    
            % circular variance
            c = sum(cos(angle(z)));
            s = sum(sin(angle(z)));
            r = sqrt(c^2 + s^2)/N;
            circVar = 1-r; 
            
        end
        
        
%         function [bgr] = local_background_single(obj,dat,varargin)
%             % LOCAL_BACKGROUND_SINGLE  under development.
%             
%             % parse input
%             pin = inputParser;
%             addOptional(pin,'method','pyfai');
%             addOptional(pin,'max_radius',10);
%             addOptional(pin,'percentile',5);
%             addOptional(pin,'bins',100);
%             addOptional(pin,'pyFAIdetector','Eiger4M');
%             addOptional(pin,'pyFAImask',[]);
%             parse(pin,varargin{:});
%             
%              % Process mask
%             mask = obj.process(obj.mask);
%             mask = round(mask ./ (obj.binx*obj.biny));
%             
%             tic;
%             switch pin.Results.method
%                 case 'zaefferer'
%                     zaeff.bgr = 0;
%                 case 'pyfai'
%                     % import modules
%                     py.importlib.import_module('numpy');
%                     py.importlib.import_module('pyFAI');
% 
%                     dimensions = py.tuple({obj.Ny,obj.Nz});
%                     bins = pin.Results.bins;
% 
%                     kwa1 = pyargs('detector',pin.Results.pyFAIdetector,'pixel1',obj.pixelsize,'pixel2',obj.pixelsize);
% 
%                     % apply special pyfai mask
%                     if ~isempty(pin.Results.pyFAImask)
%                         % reshape to fit
%                         signal1d = reshape(pin.Results.pyFAImask,1,[]);
%                         reshaped = py.numpy.reshape(signal1d,dimensions);
%                         kwa2 = pyargs('correctSolidAngle',0,'method','numpy','unit','r_mm','error_model','poisson','mask',reshaped);
%                     else
%                         kwa2 = pyargs('correctSolidAngle',0,'method','numpy','unit','r_mm','error_model','poisson');
%                     end
% 
%                     % set azimuthal integrator object
%                     ai = py.pyFAI.AzimuthalIntegrator(kwa1);
%                     ai.setFit2D(1,obj.pbz,obj.pby);
%                     
%                     % reshape to fit
%                     signal1d = reshape(dat,1,[]);
%                     reshaped = py.numpy.reshape(signal1d,dimensions);
% 
%                     tmp = ai.integrate2d(reshaped, bins, kwa2);
%                     tmp2 = py.numpy.reshape(tmp{1},py.tuple({360*bins,}));
%                     dat_2d = reshape(cell2mat(cell(tmp2.tolist())),bins,360);
%                     mm_scale = obj.q_of_n(cell2mat(cell(tmp{2}.tolist())));
%                     
%                     % apply percentile criterium
%                     filtered = dat_2d;
%                     filteredAverage = zeros(size(dat_2d,1),1);
%                     for ii = 1:size(dat_2d,1)
%                         
%                         line = dat_2d(ii,:);
%                         upper = prctile(line,100 - pin.Results.percentile);
%                         lower = prctile(line,pin.Results.percentile);
%                         fprintf(1,'lower: %f, upper: %f, percentile: %d\n',lower,upper,pin.Results.percentile);
% 
%                         line(line < lower) = nan;
%                         line(line > upper) = nan; 
%                         
%                         filtered(ii,:) = line;
%                         filteredAverage(ii) = mean(line,'omitnan');
%                     end
%                     
%                     % create artificial background by interpolating back
%                     % onto 2d detector
%                     R_mm = obj.R.*obj.pixelsize*1000;
%                     
%                     bgr = zeros(Nz,Ny);
%                     for ii = 1:obj.Nz
%                         for jj = 1:obj.Ny
%                             % find the according bin
%                             [~,I] = obj.helper.getIndex(mm_scale,R_mm(ii,jj));
%                             bgr(ii,jj) = filteredAverage(I);
%                         end
%                     end
%                     
%                     pyfai_res.bgr = bgr;
%                     pyfai_res.dat_2d = dat_2d;
%                     pyfai_res.filtered = filtered;
%                     pyfai_res.filteredAverage = filteredAverage;
%             end
%             toc;
%             
%             % return data
%             switch pin.Results.method
%                 case 'zaefferer'
%                     bgr = zaeff;
%                 case 'pyfai'
%                     bgr = pyfai_res;
%             end
%             bgr.mask = mask;
%         end
        
        
%         function [result] = remove_background_saxs(~,aavgResult1,aavgResult2,scalingFactor)
%             % REMOVE_BACKGROUND_SAXS  under development.
%             
%             aavgResult1.dat_1d = aavgResult1.dat_1d - aavgResult2.dat_1d;
%             
%             if scalingFactor ~= 0
%                 aavgResult1.dat_1d = aavgResult1.dat_1d.*(aavgResult1.qr.^scalingFactor);
%             end
%             
%             result = aavgResult1;
%         end
        
        
%         function peak_picker(obj,fh,data,varargin)
%             % PEAK_PICKER  under development.
%             %
%             %   PEAK_PICKER(col,row) finds peak.
%             %
%             %   The following parameters are accepted:
%             %
%             %       fh:: [] 
%             %          row in scan (z axis).
%             %
%             %       data:: [] 
%             %          col in scan (y axis).
%             
%             % parse input
%             pin = inputParser;
%             addOptional(pin,'sny',size(data,2));
%             addOptional(pin,'snz',size(data,1));
%             addOptional(pin,'yCrop',[1 size(data,2)]);
%             addOptional(pin,'zCrop',[1 size(data,1)]);
%             addOptional(pin,'ySkip',1);
%             addOptional(pin,'zSkip',1);
%             addOptional(pin,'yLim',0);
%             addOptional(pin,'zLim',0);
%             addOptional(pin,'emptySub','off');
%             addOptional(pin,'columnFirst','off');
%             addOptional(pin,'beamline','id13'); % not used at all
%             addOptional(pin,'scanNo',1); % so far, only stxm should be used
%             addOptional(pin,'type','stxm'); % so far, only stxm should be used
%             addOptional(pin,'orientation',0); % so far, only stxm should be used
%             parse(pin,varargin{:});
%             
%             xaxis = pin.Results.yCrop(1):pin.Results.ySkip:pin.Results.yCrop(2);
%             yaxis = pin.Results.zCrop(1):pin.Results.zSkip:pin.Results.zCrop(2);
%             
%             [x,y] = ginput(1);
%             x = round(x); y = round(y);
%             disp(['You clicked on x: ' num2str(x) '(' num2str(size(data,2)) ')' ' y: ' num2str(y) '(' num2str(size(data,1)) ')']);
%             
%             if x < 1 || y < 1 || x > size(data,2) || y > size(data,1)
%                 return
%             end
%             
%             % find the maximum in a 25 x 25 area around this position
%             subregion = data((y-12):(y+12),(x-12):(x+12));
% %             [~,max_pos] = max(obj.helper.to1dcol());
% %             [cm,cm2] = CM_nano(data((x-12):(x+12),(y-12):(y+12)));
%             [v,ind]=max(subregion);
%             [v1,ind1]=max(max(subregion));
%             Y = ind(ind1);
%             X = ind1;
% 
%             cen_y = x+X-13;
%             cen_x = y+Y-13;
%             phi = atan2(cen_y - obj.pbz, cen_x - obj.pby)*180/pi;
%             if phi < 0
%                 phi = phi + 360;
%             end
%             disp(['Center at : ' num2str(cen_x) ' and ' num2str(cen_y)]);
%             disp(['            ' num2str(obj.qyaxis(round(cen_x))) ' and ' num2str(obj.qzaxis(round(cen_y)))]);
%             disp(['Angle of rotation: ' num2str(phi)]);
%             disp(['qr position: ' num2str(obj.qr(cen_y,cen_x))]);
%                 
%             % draw circle around position
%             phi = [0:360]*pi/180;
%             r = 12;
%             xCirc = r*cos(phi);yCirc = r*sin(phi);
% 
%             figure(fh);
%             hold on;
%             line(yCirc+cen_y,...
%                  xCirc+cen_x,'color','red','LineWidth',2);
%             hold off;
%             drawnow;
%             
%             obj.peak_picker(fh,data,varargin);
%         end
        
        function ind = sub2ind(obj,y,z)
            % SUB2IND  calculates the linear index within the scan, based
            % on the row and column index.
            %
            %   LININDEX = SUB2IND(col,row) calculates linear index.
            %
            %   The following parameters are accepted:
            %
            %       row:: [] 
            %          row in scan (z axis).
            %
            %       col:: [] 
            %          col in scan (y axis).
            roi = zeros(obj.scan.SNz,obj.scan.SNy);
            roi(z,y) = 1;
            indexlist = find(flipud(rot90(roi)));
            ind = indexlist(1);
        end
        
        
        function [parameter] = moment_from_selection(obj,varargin)
            
            if nargin > 1
                method = varargin{1};
            else 
                method = 'mean';
            end
            
            % get current figure
            figure(gcf);
            axhandle = get(gca,'Children');
            data = axhandle.CData;
            
            % seelct roi
            sel = roipoly;
            data = obj.helper.to1dcol(data(sel));
            
            switch method
                case 'mean'
                    parameter = mean(data);
                case 'median'
                    parameter = median(data);
                case 'std'
                    parameter = std(data);
                case 'var'
                    parameter = var(data);

                otherwise
                    error(printf('Please provide a valid method. \n See help nanodiffraction.moment_from_selection for more information'));
            end
            
            disp(parameter);
            
        end
        
        function sub = ind2sub(obj,ind)
            % IND2SUB  calculates the row and column index based on the
            % linear index in the current scan.
            %
            %   [ROW,COL] = READ(LININDEX) calculates row and column index.
            %
            %   The following parameters are accepted:
            %
            %       linindex:: [] 
            %          Linear index (first frame: 1).
            %
            row = floor(double(ind-1)/double(obj.scan.SNy))+1; % 1 -> SNz
            col = mod((ind-1),double(obj.scan.SNy))+1; % 1 -> SNy
            sub = [row,col];
        end
        
        function clicktool(obj,vis,map,varargin)
            % CLICKTOOL  Once a parameter map of the data is shown in the
            % active figure, clicktool can be called and diffraction
            % pattern or 1d saxs curves are loaded, based on the location
            % that was clicked.
            %
            %   CLICKTOOL(vis,map,optional) reads single frame.
            %
            %   The following parameters are accepted:
            %
            %       vis:: [] 
            %          visualization tool (required).
            %
            %       map:: [] 
            %          Depending on the crop and skip parameters used
            %          during the last scan analysis, the click location
            %          might have to be mapped onto the corresponding
            %          position in the scan. It map is empty, 1x1 mapping
            %          will be used, otherwise, map has to be an array
            %          containing the following parameters 
            %           [yCrop zCrop ySkip zSkip]
            %
            %       optional:: [] 
            %          If no optional parameters are given, then a
            %          diffraction pattern will loaded, otherwise, an
            %          arbitrary number of structs can be given, containing
            %          an m:SNy:SNz data block (m: length of single saxs
            %          curve). The struct has to contain a data block and a
            %          scale, hence, should be passed e.g. as 
            %           struct('data',data,'scale',sc)
           
            
            % get current figure
            fig_map = gcf;
            
            % get size of data content
            [rows, cols] = size(getimage(fig_map));
            
            % read and understand input
            if ~isempty(map)
                % parse mapping
                yCrop = map(1);
                zCrop = map(2);
                ySkip = map(3);
                zSkip = map(4);
                
                map_ij = @(i,j) deal(i*zSkip - zSkip + zCrop,j*ySkip - ySkip + yCrop); % outputs actual row and column in scan
            else
                % 1 to 1 relation
                map_ij = @(i,j) deal(i,j); % outputs actual row and column in scan
            end
            
            % get data
            if nargin > 3
                data_in = varargin; % data is a cell array of structs
                
                % is it an array of structs (as you get using b1d, pyfai,
                % etc) or a 3d data block
                if isa(data_in{1},'struct')
                    
                    data_out = cell(numel(data_in),1);
                    
                    % convert all data structs to a 3d data volume
                    for ii = 1:numel(data_in)
                        data_tmp    = data_in{ii};
                        intensity   = zeros(numel(data_tmp(1,1).dat_1d),size(data_tmp,1),size(data_tmp,2));
                        qr          = zeros(numel(data_tmp(1,1).dat_1d),size(data_tmp,1),size(data_tmp,2));
                        for jj = 1:size(data_tmp,2)
                            for kk = 1:size(data_tmp,1)
                                intensity(:,kk,jj) = data_tmp(kk,jj).dat_1d;
                                qr(:,kk,jj) = data_tmp(kk,jj).qr;
                            end
                        end
                        data_out{ii} = struct('data',intensity,'scale',qr);
%                         data_in{ii}.data = intensity;
%                         data_in{ii}.scale = qr;
                    end
                end
                
                show_block = cell(numel(data_out),1);
                for ii = 1:numel(data_out)
                    show_block{ii} = struct('dat_1d',data_out{ii}.data(:,1,1),...
                                            'qr',data_out{ii}.scale(:,1,1));
                end
            else
                data_out = [];
            end
            
            % create new figure
            new_fig = figure;
            
            % activate first figure though
            figure(fig_map)
            
            while 1
                % get a single location and the type of input
                [y,z,button] = ginput(1);
                y = round(y); 
                z = round(z);
                fprintf(1,'You clicked on y: %d/%d  z: %d/%d\n',...
                    y,cols,z,rows);
                
                if y < 1 || z < 1 || y > cols || z > rows
                    break;
                end
                
                % row and col in scan
                [r,c] = map_ij(z,y);
                
                fprintf(1,'Corresponding scan coordinates: sy: %d/%d  sz: %d/%d\n',c,obj.scan.SNy,r,obj.scan.SNz);
                    
                figure(new_fig);
                try
                    switch button
                        case 1 % left mouse click (new plot)

                            clf;
                            if ~isempty(data_out)
                                for ii = 1:numel(data_out)
                                    show_block{ii} = struct('dat_1d',data_out{ii}.data(:,z,y),...
                                        'qr',data_out{ii}.scale(:,z,y));
                                end

                                % show 1d saxs plots
                                vis.saxs(show_block);
                            else
                                % show diffraction pattern
                                disp(['Reading frame #' num2str((r-1)*obj.scan.SNy + c)]);
                                vis.diffraction(obj.readm((r-1)*obj.scan.SNy + c),'process','off');
                            end

                        case 97 % keyboard button a ('add')

                            hold all
                            if ~isempty(data_out)
                                for ii = 1:numel(data_out)
                                    show_block{ii} = struct('dat_1d',data_out{ii}.data(:,z,y),...
                                        'qr',data_out{ii}.scale(:,z,y));
                                end

                                % show 1d saxs plots
                                vis.saxs(show_block);
                            else
                                % no default behaviour
                            end
                            hold off
                    end
                catch e
                    rethrow(e);
                end
                
                % switch back to old figure
                figure(fig_map);
            end
        end
        
        
        function data = read(obj,fn)
            % READ  reads data from data module. Note, that data is not
            % treated in any way (binning, cropping, etc.)
            %
            %   DATA = READ(fn) reads single frame.
            %
            %   The following parameters are accepted:
            %
            %       fn:: [] 
            %          Frame number (>0).
            %
            data = obj.data.read(fn);
        end
        
        function data = readp(obj,fn)
            % READP  reads data from data module. Note, that data is
            % processed using the built-in function "process".
            %
            %   DATA = READP(fn) reads single frame and processes data.
            %
            %   The following parameters are accepted:
            %
            %       fn:: [] 
            %          Frame number (>0).
            %
            data = obj.process(obj.data.read(fn));
        end
        
        function data = readm(obj,fn)
            % READM  reads data from data module. Note, that data is
            % processed using the built-in function "process" and the
            % detector mask is applied on the data frame.
            %
            %   DATA = READM(fn) reads single frame, applies detector mask
            %   and processes data.
            %
            %   The following parameters are accepted:
            %
            %       fn:: [] 
            %          Frame number (>0).
            %
            
            % Process mask
            if min(size(obj.mask) ~= [obj.Nz,obj.Ny])
                try 
                    obj.mask = obj.process_mask(obj.mask);
                    mask = 1 - obj.mask;
                catch e
                    error('Detector mask size does not match detector size');
                end
            else
                mask = 1 - obj.mask;
            end
            
            data = obj.process(obj.data.read(fn)).*mask;
        end
        
        
        function [symm_mask, mask_orig, mask, good_zero, bad_zero] = make_mask_symmetric(obj,mask)
            % MAKE_MASK_SYMMETRIC  reads in a mask and then rotates the
            % mask 180 degrees around the primary beam position. The 
            % rotated mask is compared to the unrotated mask and logically
            % combined (&).
            %
            %   SYMM_MASK = MAKE_MASK_SYMMETRIC(MASK)
            %
            %   The following parameters are accepted:
            %
            %       mask:: [] 
            %          A valid detector mask.
            %            
            Ny = obj.Ny;
            Nz = obj.Nz;
            pby = round(obj.pby);
            pbz = round(obj.pbz);
            
            % rotate 180 degrees
            mask_rot = rot90(mask,2);
            
            if pby < Ny/2 && pbz < Nz/2
                % top-left
                y1 = 2*pby - 1;
                z1 = 2*pbz - 1;
                y2 = Ny - y1;
                z2 = Nz - z1;
            elseif pby > Ny/2 && pbz < Nz/2
                % top-right
                y2 = 2*(Ny - pby) + 1;
                y1 = Ny - y2;
                z1 = 2*pbz - 1;
                z2 = Nz - z1;
            elseif pby < Ny/2 && pbz > Nz/2
                % bottom-left
                y1 = 2*pby - 1;
                z2 = 2*(Nz - pbz) + 1;
                y2 = Ny - y1;
                z1 = Nz - z2;
            else
                % bottom-right
                y2 = 2*(Ny - pby) + 1;
                y1 = Ny - y2;
                z2 = 2*(Nz - pbz) + 1;
                z1 = Nz - z2;
            end
            
            
            % case comparison
            mask_orig = mask;
            bad_zero = zeros(size(mask));
            good_zero = zeros(size(mask));
            if pby < Ny/2 && pbz < Nz/2
                % upper-left
                mask_orig(1:z1,1:y1) = mask_rot(z2+1:end,y2+1:end);
                bad_zero(1:z1,1:y1) = mask(1:z1,1:y1);
                good_zero(1:z1,1:y1) = mask_rot(z2+1:end,y2+1:end);
            elseif pby > Ny/2 && pbz < Nz/2
                % top-right
                mask_orig(1:z1,y1+1:end) = mask_rot(z2+1:end,1:y2);
                bad_zero(1:z1,y1+1:end) = mask(1:z1,y1+1:end);
                good_zero(1:z1,y1+1:end) = mask_rot(z2+1:end,1:y2);
            elseif pby < Ny/2 && pbz > Nz/2
                % bottom-left
                mask_orig(z1+1:end,1:y1) = mask_rot(1:z2,y2+1:end);
                bad_zero(z1+1:end,1:y1) = mask(z1+1:end,1:y1);
                good_zero(z1+1:end,1:y1) = mask_rot(1:z2,y2+1:end);
            else
                mask_orig(z1+1:end,y1+1:end) = mask_rot(1:z2,1:y2);
                bad_zero(z1+1:end,y1+1:end) = mask(z1+1:end,y1+1:end);
                good_zero(z1+1:end,y1+1:end) = mask_rot(1:z2,1:y2);
            end
            
            % combine masks
            symm_mask = mask | mask_orig;

        end

        function [image] = azimuthal_spread(obj,dat_1d,qr)
            % AZIMUTHAL_SPREAD  takes a one-dimensional intensity profile
            % and the corresponding q_r position as input and
            % reinterpolates the data onto the two-dimensional detector
            % grid. This can be useful when semi-transparent objects such 
            % as a beamstop shall be corrected for.
            %
            %   IMAGE = AZIMUTHAL_SPREAD(DAT_1D,QR) 
            %
            %   The following parameters are accepted:
            %
            %       dat_1d:: [] 
            %          One-dimensional structure factor as obtained e.g.
            %          from the b1d or pyfai routine
            %
            %       qr:: []
            %           Corresponding Qr positions.
            %            
            xq = obj.helper.to1dcol(obj.qr);
            vq = interp1(qr,dat_1d,xq);
            image = obj.helper.to2d(vq,obj.Ny,obj.Nz);
            
        end
        
        function stitched = stitch(obj,data,stitchy,stitchz,varargin)
            % STITCH  takes a one-dimensional cell array of two-dimensional 
            % data and combines all data into a single array, based on the
            % stitching dimensions.
            %
            %   STITCHED = STITCH(DATASET,STITCHY,STITCHZ,WINDOW) 
            %
            %   The following parameters are accepted:
            %
            %       DATASET:: [] 
            %          One-dimensional cell array data{:} that contains
            %          two-dimensional diffraction patterns of identical
            %          dimensions.
            %
            %       STITCHY:: []
            %           Number of images to be concatenated horizontally.
            %
            %       STITCHZ:: []
            %           Number of images to be concatenated vertically.
            %            
            %       WINDOW:: [] (optional)
            %           4-element vector containing the lower and upper
            %           limits (in pixel units) along the horizontal and
            %           lower and upper limits along the vertical dimension
            %           of the scattering pattern. 
            %           Example: [100 900 100 900]
            
            defaults = struct('win',[],'append',0,'prepend',0);
            fields = fieldnames(defaults);
            if nargin > 2
                opts = varargin{1};
                for f = 1:numel(fields)
                    if ~isfield(opts,fields{f})
                        opts.(fields{f}) = defaults.(fields{f});
                    end
                end
            else
                opts = defaults;
            end
            
            prepend = @(cellArr) [ {zeros(size(cellArr{1}))},cellArr];
            append = @(cellArr) [ cellArr, {zeros(size(cellArr{1}))}];
            
            if opts.append > 0
                for ii = 1:opts.append;data = append(data);end
            end
            if opts.prepend > 0
                for ii = 1:opts.prepend;data = prepend(data);end
            end
            
            if numel(data) ~= stitchy*stitchz
                [m,n] = size(data);
                warning('data matrix (%d x %d) does not match stitch size (%d x %d)\n Zeros will be added.',m,n,stitchz,stitchy);
                for ii = 1:(stitchy*stitchz - m*n);data = append(data);end
            end
            
            usewin = 0;
            if ~isempty(opts.win)
                win = opts.win;
                yl = win(1);
                yu = win(2);
                zl = win(3);
                zu = win(4);
                
                usewin = 1;
            end
                
            data_tmp = cell(stitchy*stitchz,1);
            for s = 1:stitchy*stitchz
                if usewin
                    data_tmp{s} = data{s}(zl:zu,yl:yu);
                else
                    data_tmp{s} = data{s};
                end
            end
            
            line = cell(stitchz,1);
            for l = (1:stitchz) - 1
                line{l+1} = cat(2,data_tmp{(1:stitchy) + l*stitchy});
            end
            stitched = cat(1,line{:});
        end
        
    end
end
