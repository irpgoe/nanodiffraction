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
    %   coordinates. Here is a short description of the available helper
    %   functions stored in EXPERIMENT.helper
    %    n_of_q(qr): Requires a reciprocal wavevector coordinate and
    %    calculates the according pixel coordinate
    %    q_of_n(n): Requires a pixel coordinate and calculates the
    %    according reciprocal wavevector coordinate
    %    twoTh_of_n(n): Requires a pixel coordinate and calculates the
    %    according scattering angle 2*theta
    %    to1dcol(x): Requires a 2d array x and transforms this array into a
    %    1d column vector
    %    to1drow(x): Requires a 2d array x and transforms this array into a
    %    1d row vector
    %    to2d(x,ny,nz): Requires a 1d array x and transforms this array
    %    into a 2d array with dimensions (nz,ny)
    %    reduce(i,start,skip): Simplified indexing
    %    simplify_sf(data): Not used 
    %    scale(x): Requires a 2d array x and rescales the array into the
    %    range (0,1) based on the minimum and maximum values in the array
    %    scaleH(x,h): Requires a 2d array x and rescales the array into the
    %    range (0,h) based on the minimum and maximum values in the array
    %    and the high value h
    %    scaleLH(x,l,h): Requires a 2d array x and rescales the array into
    %    the range (l,h) based on the minimum and maximum values in the
    %    array and the low value l and the high value h
    %    clipL(x,l): Requires a 1d or 2d array x. Values below the clipping
    %    value l will be set to l.
    %    clipH(x,h): Requires a 1d or 2d array x. Values above the clipping
    %    value h will be set to h.
    %    clipLH(x,l,h): Requires a 1d or 2d array x. Values below the 
    %    clipping value l will be set to l and values above the clipping
    %    value h will be set to h.
    %    autoClip(x): Requires a 1d or 2d array and calculates the 5%
    %    low and high percentile. 
    %    getIndex(x,xloc): Calculates the closest linear index in the 1d 
    %    array x based on the query value xloc
    %    addTransp(transp): Overlays a transparency map over the current
    %    figure. Values in the transparency map should range between (0,1)
    %    and the size of the transparency map has to match the size of the
    %    data in the figure
    %    normDist(x): Removes the mean and divides by the standard
    %    deviation of a distribution x.
    %
    %   In order to perform calculations on real data, one has to couple
    %   the nanodiffraction class to the files class to "tell the
    %   experiment, where the data is stored". This can be achieved by
    %   setting EXPERIMENT.data = [handle to files class] or by using the
    %   helper function link([handle to files class],[handle to helper
    %   class]). See help link for more information.
    %
    %   Furthermore, three other properties of the nanodiffraction class
    %   are of interest:
    %       - nanodiffraction.mask:
    %           A bad pixel mask can be set using the function set_mask().
    %       - nanodiffraction.sels
    %           A mask identifying pixels that shall be analyzed can be
    %           defined using set_selection(). More selections can be added
    %           using add_selection(). This is useful since data has to be
    %           read only once, while the analysis can be performed
    %           multiple times on the same data (but different selections).
    %           In many cases it is useful to use the functions
    %           radial_mask() and azimuthal_mask() in order to define
    %           selections. 
    %       - nanodiffraction.corr
    %           A correction matrix that will be multiplied with the data
    %           prior to the analysis to compensate for semi-transparent
    %           objects such as glass capillaries.
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
    %     files-class handle. This handle can also be set using the link()
    %     function. 
    %
    %   Example::
    %      experiment =
    %      nanodiffraction('energy',12.8e3,'detDistance',1.2,...
    %      'pby',121.2,'pby',144.1,...
    %      'Ny',1100,'Nz',1220,...
    %      'fluoCounters',{{'total',1 4100},{'caka',600,660}});
    %
    %   Methods::
    %      See seperate help for each method. Also, refer to the methods
    %      summary in >> doc nanodiffraction.
    
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
        data = files;       % holds link to the data
        helper = struct;    % holds helper functions
        scan = struct;      % holds scan parameters
        
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
        phi2 = [];
        
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
            addOptional(pin,'energy',13.8e3);
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
            obj.phi2 = atan2d(repmat(obj.zaxis',1,obj.Ny),repmat(obj.yaxis,obj.Nz,1));
            obj.phi2(obj.phi2<0) = obj.phi2(obj.phi2<0) + 180;
            
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
        
        
        function p = get_parameters(obj)
            % GET_PARAMETERS  Returns the experimental parameters stored in
            % the nanodiffraction class. Note, that helper functions and
            % links to other classes are not returned.
            %
            %   result = get_parameters()
            %
            % The following arguments are supported:
            %   This function does not require any arguments.
            %
            % Example:
            %   e = nanodiffraction();
            %   p = e.get_parameters()
            %
            % Output arguments:
            %   result:: Structure containing all experimental parameters
            %       stored in this instance of the nanodiffraction class.
            %
            
            p = struct();
            pList = properties(obj);
            for i = 3:numel(pList)
                p.(pList{i}) = obj.(pList{i});
            end
        end
        
        
        function result = get_location_in_frame(obj,y,z,sn,stitch)
            % GET_LOCATION_IN_FRAME  Calculates the scan position within a
            % frame based on the location in the stitched image.
            %
            %   result = get_location_in_frame(y,z,sn,stitch);
            %
            % The following arguments are supported:
            %   y:: [] (required)
            %       Description.
            %
            %   z:: [] (required)
            %       Description.
            %
            %   sn:: [] (required)
            %       Description.
            %
            %   stitch:: [] (required)
            %       Description.
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %	result:: Structure containing the following fields
            %       - linind:: Description.
            %       - linfield:: Description.
            %       - field:: Description.
            %       - ind:: Description.
            %

            sny = sn(1);
            snz = sn(2);
            stitchy = stitch(1);
            stitchz = stitch(2);
            
            if size(y,1) > 1
                y = y';
            end
            if size(z,1) > 1
                z = z';
            end
            
            result.ind = [y-floor((y-1)/sny)*sny, z-floor((z-1)/snz)*snz];
            result.linind = (result.ind(2)-1)*sny + result.ind(1);
            result.field = [floor((y-1)/sny) + 1, floor((z-1)/snz) + 1];
            result.linfield = result.field(1) + (result.field(2)-1)*stitchy;
           
        end
        
        
        function [varargout] = get_index(obj,q,pos)
            % GET_INDEX  Calculates pixel coordinate in pixels given a Q
            % value.
            %
            %   out = get_index(q,position)
            %
            % The following arguments are supported:
            %   q:: [] (required)
            %       q-axis, where q is the reciprocal wavevector coordinate
            %       in units of reciprocal nanometers.
            %
            %   position:: [] (required)
            %       Query position. 
            %
            % Example: 
            %   Example missing.
            %
            % Output arguments:
            %   out:: Linear index of the nearest matching element in the
            %   q-axis.
            %
            
            for i = 1:numel(pos)
                [~,varargout{i}] = obj.helper.getIndex(q,pos(i));
            end
        end
        
        function [n] = n_of_q(obj,q)
            % N_OF_Q  Calculates pixel coordinate in pixels given a Q
            % value. This function cannot be defined as an anonymous 
            % function, because pixelsize might in fact not be constant 
            % (and the other constants in some rare occasions as well)
            %
            %   n = n_of_q(q)
            %
            % The following arguments are supported:
            %   q:: [] (required)
            %       Coordinate in units of the reciprocal wavevector
            %       (inverse nanometers)
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   n:: Coordinate in units of detector pixels.
            %
            
            n = tan(2*asin(q*obj.wavelength*1E9/4/pi)) * obj.detDistance/obj.pixelsize;
        end
        
        function [q] = q_of_n(obj,n)
            % Q_OF_N  Calculates Q coordinate in [nm^-1] given a pixel
            % coordinate. This function cannot be defined as an anonymous 
            % function, because pixelsize might in fact not be constant 
            % (and the other constants in some rare occasions as well)
            %
            %   q = q_of_n(n)
            %
            % The following arguments are supported:
            %   n:: [] (required)
            %       Coordinate in units of detector pixels
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   q:: Coordinate in units of the radial wavevector transfer
            %   (reciprocal nanometers).
            %
            
            q = 4*pi/obj.wavelength/1E9 * sin(atan(n*obj.pixelsize/obj.detDistance)/2);
        end
        
        
        function [q,x,y,z] = qr_tilt(obj,alpha,beta,gamma)
            % QR_TILT  Description missing.
            %
            %   [qr,x,y,z] = qr_tilt(alpha, beta, gamma) 
            %
            % The following arguments are supported:
            %   alpha:: []
            %       Tilt angle of the detector in degrees.
            %
            %   beta:: []
            %       Tilt angle of the detector in degrees.
            %
            %   gamma:: []
            %       Tilt angle of the detector in degrees.
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   qr:: Corrected radial wavevector transfer detector 
            %   coordinates, based on the detector tilt angles.
            %
            %   x:: Corrected x-position of the detector pixel
            %   coordinates, based on the detector tilt angles.
            %
            %   y:: Corrected y-position of the detector pixel
            %   coordinates, based on the detector tilt angles.
            %
            %   z:: Corrected z-position of the detector pixel
            %   coordinates, based on the detector tilt angles.
            %
            
            % define rotation matrices
            rotx = [[1 0 0];[0 cosd(alpha) -sind(alpha)];[0 sind(alpha) cosd(alpha)]];
            roty = [[cosd(beta) 0 sind(beta)];[0 1 0];[-sind(beta) 0 cosd(beta)]];
            rotz = [[cosd(gamma) -sind(gamma) 0];[sind(gamma) cosd(gamma) 0];[0 0 1]];
            
            % apply the rotations in inverse order
            rot = rotx*roty*rotz;
            
            % apply rotation on each pixel to get the correct pixel
            % coordinates
            y = repmat(obj.yaxis, obj.Nz, 1).*obj.pixelsize;
            z = repmat(transpose(obj.zaxis), 1, obj.Ny).*obj.pixelsize;
            x = repmat(0, obj.Nz, obj.Ny);
            
            q = zeros(obj.Nz, obj.Ny);
            for i = 1:obj.Nz
                for j = 1:obj.Ny
                    new_coordinates = rot * [x(i,j);y(i,j);z(i,j)];
                    x(i,j) = new_coordinates(1) + obj.detDistance;
                    y(i,j) = new_coordinates(2);
                    z(i,j) = new_coordinates(3);
            
                    % calculate qr
                    theta = atan(sqrt(y(i,j).^2 + z(i,j).^2)/x(i,j));
                    q(i,j) = 4*pi/obj.wavelength/1E9 * sin(theta/2);
                end
            end
        end
        
        
        function refresh(obj)
            % REFRESH  This function is needed in case the geometry 
            % settings change. This could happen eg. when for 
            % some reason the detector pixel size was changed (e.g. when
            % binning was used) and q_r, primary beam position, etc. has to
            % be recalculated.
            %
            %   refresh()
            %
            % The following arguments are supported:
            %   This function does not require any arguments.
            %
            % Example:
            %   e = nanodiffraction();
            %   e.pby = 1000; e.pbz = 1200; e.detDistance = 1.2;
            %   e.refresh();
            %
            % Output arguments:
            %   This function does not return any arguments.
            %
            
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
            obj.phi2 = atan2d(repmat(obj.zaxis',1,obj.Ny),repmat(obj.yaxis,obj.Nz,1));
            obj.phi2(obj.phi2<0) = obj.phi(obj.phi2<0) + 180;
            
            % detector roi in integers
            obj.detRoiY = uint32(obj.detRoiY);
            obj.detRoiZ = uint32(obj.detRoiZ);
            
            % update masks
            if max(size(obj.mask) ~= [obj.Nz,obj.Ny])
                fprintf(1,'Processing detector mask...\n');
                try 
                    obj.mask = obj.process_mask(obj.mask);
                catch e
                    try 
                        obj.mask = obj.process_mask(obj.mask_orig);
                    catch err
                        error(err);
                    end
                end
            end
            if max(size(obj.sels(:,:,1)) ~= [obj.Nz,obj.Ny])
                fprintf(1,'Processing selections...\n');
                tmp = zeros(obj.Nz,obj.Ny,size(obj.sels,3));
                for ii = 1:size(obj.sels,3)
                    try 
                        tmp(:,:,ii) = obj.process_mask(obj.sels(:,:,ii));
                    catch err
                        try 
                            tmp(:,:,ii) = obj.process_mask(obj.sels_orig(:,:,ii));
                        catch err
                            error(err);
                        end
                    end
%                     tmp(:,:,ii) = obj.process_mask(obj.sels(:,:,ii));
                end
                obj.sels = tmp;
            end     
            if max(size(obj.corr) ~= [obj.Nz,obj.Ny])
                fprintf(1,'Processing correction matrix...\n');
                try 
                    obj.corr = obj.process_mask(obj.corr);
                catch err
                    try 
                        obj.corr = obj.process_mask(obj.corr_orig);
                    catch err
                        error(err);
                    end
                end
            end
        end
        
        
        function new = copy(obj)
            % COPY  Makes a copy of a handle object. This keeps track of 
            % all the given parameters so that multiple instances of a 
            % display module can be used.
            %
            %   new_instance = copy(nanodiffraction_object)
            %
            % The following arguments are supported:
            %   
            %   nanodiffraction_object:: [] (required)
            %       This function requires a nanodiffraction object that 
            %       should be copied.
            %
            % Example:
            %   e_generic = nanodiffraction();
            %   e_saxsdata = copy(e_generic);
            %   e_waxsdata = copy(e_generic);
            %
            % Output arguments:
            %
            %   This function outputs the following arguments:
            %       
            %       new_instance:: A copy of the initial 
            %       nanodiffraction_object.
            %
            
            
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
            % be chosen in a random order and combination. However,
            % corrections are always performed before any subsequent 
            % analysis. The "healing" functionality has a special role, as
            % it is a method but is always performed before any other
            % analysis but after background and other corrections, for
            % obvious reasons.
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
            %        indicates a scan point that will be analyzed.
            %
            %      yCrop, zCrop:: [start end] 
            %        Start and end can be given to define a certain
            %        range within the scan to be analyzed.
            % 
            %      ySkip, zSkip:: [1]
            %        If only every n-th scan point should be analyzed
            %
            %      emptySub:: ['off'] 
            %        Subtracts an empty image from the data immediately 
            %        after the data has been read. If an empty image is
            %        given (see keyword 'empty'), it will be used,
            %        otherwise, it will be attempted to use the empty image
            %        that was defined at an earlier stage using
            %        set_empty().
            %
            %        Options: 'on' | 'off'
            %
            %      correction:: ['off'] 
            %        Multiplies a 2d
            %        correction matrix with the data. E.g. used to account
            %        for semitransparent objects obscuring the detector
            %
            %        Options: 'on' | 'off'
            %
            %      method:: ['none'] 
            %        A string composed of the following options: 'debug', 
            %        'rotate','pca','stxm','maxProj','crystal','fluo',
            %        'average','sum','data','circDist'. 
            %        Methods can be combined using a '+' symbol between the
            %        methods:
            %           e.g. 'stxm+pca+maxProj+sum'
            %
            %      fluoCounters:: [{}] 
            %         The format is {{<name>,<channel_low>,<channel_high>},...}, 
            %           e.g. {{'caka',1000,1400},{},...}.
            %
            %      parallel:: ['off']
            %        When turned on, all available workers will be used to
            %        speed up the calculation process. 
            %        Warning: This is currently only implemented for the 
            %        stxm method.
            %
            %      scanmode:: ['horz']
            %        See description below.
            %
            %        Options: 'horz'|'horzalt'|'vert'|'vertalt'
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
            %           - 'qy':   2d array of horizontal wavevector transfer [repmat(obj.qyaxis ,obj.Nz,1)]
            %           - 'qz':   2d array of vertical wavevector transfer [repmat(obj.qzaxis',1,obj.Ny)]
            %           - 'sel':  Logical 2d array. A value of 1 represents a pixel that will be analyzed [obj.sels]
            %           - 'mask': Logical 2d array. A value of 1 represents a pixel that will be discarded [obj.mask]
            %
            %   stxm: Integrates the scattered intensity of each
            %   diffraction pattern (darkfield contrast).
            %
            %       PARAMETER ARGUMENTS: 'stxmArgs',struct()
            %           struct() can take the following arguments:
            %           - 'qrMin':  Only pixels above the radial wavevector transfer qrMin will be analyzed [-1]
            %           - 'qrMax':  Only pixels below the radial wavevector transfer qrMax will be analyzed [-1]
            %           - 'sel':    Logical 2d array. A value of 1 represents a pixel that will be analyzed [ones(obj.Nz,oabj.Ny))]
            %           - 'mask':   Logical 2d array. A value of 1 represents a pixel that will be discarded [obj.mask]                     
            %            
            %   circDist: Calculates the circular mean and circular'grid',obj.phi2,'bins',360,'sel',obj.sels,'mask',obj.mask,'bgr',[]),...
            %   variance of a distribution. Here, used to calculate the
            %   circular mean of I(phi) for a given qr-Interval.
            %
            %       PARAMETER ARGUMENTS: 'circDistArgs',struct()
            %           struct() can take the following arguments:
            %           - 'grid': Pixel map of azimuthal coordinates [obj.phi2]
            %           - 'bins': Number of azimuthal bins [360]
            %           - 'bgr':  2d background image that is going to be subtracted before the mean scattering angle is determined []
            %           - 'sel':  Logical 2d array. A value of 1 represents a pixel that will be analyzed [obj.sels]
            %           - 'mask': Logical 2d array. A value of 1 represents a pixel that will be discarded [obj.mask]                        
            %
            %   maxProjection: Calculates the maximum projection of a stack of
            %   diffraction patterns. No arguments are required.
            %
            %   maxProj: Synonymous to maxProjection
            %
            %   mp: Synonymous to maxProjection
            %
            %   crystal: Performs a thresholding for crystal mapping based
            %   on SNR criterium.
            %
            %       PARAMETER ARGUMENTS: 'crystalArgs',struct()
            %           struct() can take the following arguments:
            %           - 'bgr':       2d background image that is going to be subtracted before the summing procedure [zeros(obj.Nz,obj.Ny)]
            %           - 'threshold': All pixels with counts above threshold will be summed [50]
            %           - 'qrMin':     Only pixels above the radial wavevector transfer qrMin will be analyzed [0]
            %           - 'qrMax':     Only pixels below the radial wavevector transfer qrMax will be analyzed [100]
            %           - 'sel':       Logical 2d array. A value of 1 represents a pixel that will be analyzed [obj.sels]
            %           - 'mask':      Logical 2d array. A value of 1 represents a pixel that will be discarded [obj.mask]                     
            %
            %   pyfai: Fast azimuthal integration based on the python pyfai
            %   library.
            %
            %       PARAMETER ARGUMENTS: 'pyfaiArgs',struct()
            %           struct() can take the following arguments:
            %           - 'bins': Number of bins along radial axis [100]
            %           - 'mask': Logical 2d array. A value of 1 represents a pixel that will be discarded []
            %           - 'cake': Toggle cake integration, can be either 'on' or 'off' [off]
            %           - 'cakeRoi': If 'cake' is 'on', then cakeRoi can be set within a range (phi1,phi2), where phi1 and phi2 are angles in degrees within the range (0,360) and phi2 should be smaller than phi1 [0 360]            
            %           - 'angles': []
            %           - 'to2D': Can be either 'on' or 'off'. If 'on', then the 2d regridded pattern is output ['off']
            %           - 'roi2d': ['off']
            %
            %   b1d: Fast azimuthal integration based on rebinning the
            %   data.
            %
            %       PARAMETER ARGUMENTS: 'b1dArgs',struct()
            %           struct() can take the following arguments:
            %           - 'grid': Pixel coordinates, e.g. a 2d map of the radial wavevector transfer [obj.qr]
            %           - 'bins': Number of bins along radial axis [100]
            %           - 'sel':  Logical 2d array. A value of 1 represents a pixel that will be analyzed [ones(obj.Nz,obj.Ny)]
            %           - 'mask': Logical 2d array. A value of 1 represents a pixel that will be discarded [obj.mask]
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
            %           - 'bgr': 2d background image that is going to be subtracted before the summing procedure [zeros(obj.Nz,obj.Ny)]
            %           - 'sigma': Parameter of Gaussian filter [1.6]
            %           - 'noise_level': SNR level that functions as a threshold for peak detection [6]
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
            %        *  v  v  v  v
            %        v  v  v  v  .
            %        v  v  v  v  .
            %        v  v  v  v  .
            %     4 vertical,alternating <vertalt>
            %        *   A > v   A > v
            %        v   A   v   A   .
            %        v   A   v   A   .
            %        v > A   v > A   .
            %
            %   Output arguments:
            %       Result:: 
            %           struct containing the results for each method.
            %           Syntax: struct_name.method_name.method_results
            %
            %       Current options:
            %           struct_name.
            %               data (processed data stack)
            %               maxProjection (maximum projection)
            %               mp (same as maxProjection)
            %               crystal (crystal map)
            %               devel
            %               dpc - dpcx (differential phase contrast,
            %                     horizontal)
            %                     dpcy (differential phase contrast,
            %                     vertical)
            %               stxm - df (darkfield)
            %               pca - w (anisotropy)
            %                     angle (orientation of the scattering)
            %               fluo
            %               avg
            %               sum
            %               rotate
            %               symmetry - Ln
            %               b1d - dat_1d (azimuthally averaged intensities)
            %                     qr (radial axis)
            %                     error (error metric, L2 norm)
            %               pyfai - dat_1d (azimuthally averaged intensities)
            %                       r (radial axis)
            %                       error (error metric, L2 norm)
            %
            %       Note, that the output structure can be split into
            %       variables using the helper function
            %       split_struct(). This should improve readability
            %       of the code.
            
            
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
            addOptional(pin,'lineEmpty',[]);
            addOptional(pin,'correction','off',@ischar);
            addOptional(pin,'parallel','off',@ischar);
            addOptional(pin,'live','off',@ischar);
            % stxm, crystal, b1d and symmetry
            addOptional(pin,'liveArgs',struct(),@isstruct);
            addOptional(pin,'stxmArgs',struct(),@isstruct);
            addOptional(pin,'crystalArgs',struct(),@isstruct);
            addOptional(pin,'symmetryArgs',struct(),@isstruct);
            addOptional(pin,'circDistArgs',struct(),@isstruct);
            addOptional(pin,'b1dArgs',struct(),@isstruct);
            addOptional(pin,'pcaArgs',struct(),@isstruct);
            addOptional(pin,'pyfaiArgs',struct(),@isstruct);
            addOptional(pin,'avgPerLineArgs',struct(),@isstruct);
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
            lineEmptyYes = 0;
            if strcmp(pin.Results.emptySub,'on')
                if ~isempty(pin.Results.empty)
                    emptyYes = 1;
                    empty = pin.Results.empty;
                elseif ~isempty(pin.Results.lineEmpty)
                    lineEmptyYes = 1;
                    lineEmpty = pin.Results.lineEmpty; 
                    empty = lineEmpty{1};   
                else
                    emptyYes = 1;
                    empty = obj.empty;
                end
            else
                empty = obj.empty;
            end

            % verify all masks and selections
            mask    = obj.verify_array(obj.mask, logical(obj.mask), @obj.process_mask);
            corr    = obj.verify_array(obj.corr, obj.corr,          @obj.process);
            sels    = obj.verify_array(obj.sels, obj.sels,          @obj.process);
            empty   = obj.verify_array(empty,    empty,             @obj.process);
                
            roiYes = 0;
            if ~isempty(pin.Results.roi)
                roiYes = 1;
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
            itotYes = 0;
            imaxYes = 0;
            circDistYes = 0;
            avgPerLineYes = 0;
            
            % Slice some variables
            fnrs = obj.scan.fnrs;
            
            
            % All methods are initialized with default values if they are
            % selected using the 'method' parameter
            nr_rows = length(zsampling);
            nr_cols = length(ysampling);
            if contains(method,'circDist')
                circDistYes = 1;
                phi = zeros(nr_rows*nr_cols,1);
                w = zeros(nr_rows*nr_cols,1);
                circDistArgs = update_defaults(...
                    struct('grid',obj.phi2,'bins',360,'sel',obj.sels,'mask',obj.mask,'bgr',[]),...
                    pin.Results.circDistArgs);
                
                [ph,sortIndex_rad,n_rad,bin_rad] = prepare_averaging(circDistArgs.grid,circDistArgs.mask,circDistArgs.sel,circDistArgs.bins);
            end
            if contains(method,'itot')
                itotYes = 1;
                itot = zeros(nr_rows*nr_cols,1);
            end
            if contains(method,'imax')
                imaxYes = 1;
                imax = zeros(nr_rows*nr_cols,1);
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
                phi = zeros(nr_rows*nr_cols,size(sels,3));
                w = zeros(nr_rows*nr_cols,size(sels,3));
                pcaArgs = update_defaults(...
                    struct('qy',repmat(obj.qyaxis ,obj.Nz,1),'qz',repmat(obj.qzaxis',1,obj.Ny),'sel',sels,'mask',mask),...
                    pin.Results.pcaArgs);
            end
            if contains(method,'stxm') 
                stxmYes = 1;
                stxmArgs = update_defaults(...
                    struct('mask',obj.mask,'qrMin',-1,'qrMax',-1,'sel',sels),...
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
            if contains(method,'maxProj') || contains(method,'mp') || contains(method,'maxProjection')% legacy
                maxProjYes = 1;
                maxproj = zeros(obj.Nz,obj.Ny);
            end
            if contains(method,'crystal')
                crystalYes = 1;
                crystal = zeros(nr_rows*nr_cols,1);
                crystalArgs = update_defaults(...
                    struct('bgr',zeros(obj.Nz,obj.Ny),'mask',mask,'threshold',50,'qrMin',0,'qrMax',100,'sel',sels),...
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
            if contains(method,'average') || contains(method,'avg')
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
                Ln = zeros(nr_rows*nr_cols,length(0:4:ceil(max(obj.R(:)))));               
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
                    struct('grid',obj.qr,'bins',100,'mask',mask,'sel',sels,'bgr',[]),...
                    pin.Results.b1dArgs);
                
                % prepare averaging requires only a single mask and
                % selection
                if size(b1dArgs.sel,3) > 1 
                    warning('Only slower method can be used at the moment')
                else
                    [qr,sortIndex_azim,n_azim,bin_azim] = prepare_averaging(b1dArgs.grid,b1dArgs.mask,b1dArgs.sel,b1dArgs.bins);
                end
            end
            if contains(method,'fluo') 
                fluoYes = 1;
                % initialize struct fields for results corresponding to
                fluoResult = struct;
                fluoResult.averageSpectrum = zeros(4096,1);
                fluoResult.sumSpectrum = zeros(4096,1);
                for ii = 1:length(obj.fluoCounters)
                    tmp = obj.fluoCounters{ii}(1);
                    % register field
                    fluoResult.(tmp{1}) = zeros(nr_rows * nr_cols,1);
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
            
            % after the initialization, live plotting
            if strcmp(pin.Results.live,'on')
                fprintf(1,'Information on live plotting: The first selected method will be shown live');
                liveYes = 1;
                line = 0;
                liveMethod = strsplit(method,'+');
                liveMethod = liveMethod{1};
                methodKeys = {'stxm','pca','circDist','maxProj','average','Itot','Imax','crystal'};
                methodVarNames = {'df','phi','phi','maxproj','sumResult','itot','imax','crystal'};
                methodMap = containers.Map(methodKeys,methodVarNames);
                try 
                    methodVariable = methodMap(liveMethod);
                    fprintf(1,'You have chosen the following method: %s\n',liveMethod);
                    eval(['methodVariable = ' methodVariable ';']);
                catch err
                    warning('Method not supported for live display, live mode will be disabled')
                    liveYes = 0;
                end
            end
            if contains(method,'live')
                liveYes = 1;
                line = 0;
                displayClass = display();
                link(obj,displayClass);
                liveArgs = update_defaults(...
                    struct('display',displayClass,'cAxis',[0 1]),...
                    pin.Results.liveArgs);
            end
            if contains(method,'avgPerLine')
                avgPerLineYes = 1;
                
                % initialize cell array
                avgPerLine = cell(obj.scan.SNz,1);
                [avgPerLine{1:obj.scan.SNz, 1}] = deal(zeros(obj.Nz,obj.Ny));
                                
                avgPerLineArgs = update_defaults(...
                    struct('maxFramesPerLine',100),...
                    pin.Results.avgPerLineArgs);
                    
                roi_buff= zeros(size(roi));
                for l = 1:size(roi,1)
                    tmp = find(roi(l,:));
                    if length(tmp) == 0
                         warning(['Line ' num2str(l) ' is empty.']);
                    elseif length(tmp) >= avgPerLineArgs.maxFramesPerLine
                        roi_buff(l,tmp(1:avgPerLineArgs.maxFramesPerLine)) = 1;
                    else
                        roi_buff(l,tmp(1:end)) = 1;
                    end
                end
                roi = roi_buff;
            end
            
            
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
            
            
            tic;       
            
            if ~strcmp(pin.Results.parallel,'on')
                % serial version
                fprintf(1,'%6.2f%%', 0.00);
                for jj = 1:length(indexlist)

                    % progress 
                    fprintf(1,'\b\b\b\b\b\b\b%6.2f%%', double(jj)/length(indexlist)*100);

                    index = indexlist(jj);  % index is column-based
                    fn = fnrs(index);       % otherwise this would not work
                    line = floor((index-1)/obj.scan.SNy) + 1;
                    
                    if roiYes 
                        ss = index;
                    else
                        ss = jj;
                    end
                    
                    % read data if necessary/requested
                    if ~contains(method,'none') && ~contains(method,'fluo')

                        dat = obj.process(obj.data.read(fn));
                        % dat = obj.read_wait(fn);
                        
                        % perform pre-processing
                        if emptyYes
                            dat = dat - empty;
                        end

                        if lineEmptyYes
                            dat = dat - lineEmpty{line};
%                             dat = dat - transpose(lineEmpty{line});
                        end  
                        
                        if strcmp(pin.Results.correction,'on')
                            dat = dat.*corr;
                        end
                        
                        % always apply mask
                        % dat(mask) = 0; % is this still needed?
                    end
                    
                    % only fully pre-processed images should be healed
                    if healYes
                        dat = dat(healIl);
                    end
                    
%                     % live viewing
%                     if liveYes
%                         % show line progress
%                         if index > (line+1)*numel(ysampling) % line done
%                             if strcmp(liveMethod,'maxProj') || strcmp(liveMethod,'average')
%                                 imagesc(methodVariable);axis image;colorbar;colormap(pmkmp);title(liveMethod);
%                             elseif strcmp(liveMethod,'pca') || strcmp(liveMethod,'circDist')
%                                 methodVariable_tmp = sort_scan(methodVariable,nr_rows,nr_cols,scanmode);
%                                 imagesc(methodVariable_tmp-90);colormap(cmap('C2'));colorbar;axis image;
%                             else
%                                 methodVariable_tmp = sort_scan(methodVariable,nr_rows,nr_cols,scanmode);
%                                 imagesc(methodVariable_tmp);colormap(pmkmp);colorbar;axis image;
%                             end
%                             drawnow;
%                             
%                             % wait for next line
%                             line = line + 1;
%                         end
%                         
%                     end

                    % live viewing
                    if liveYes
                        liveArgs.display.diffraction(dat.*~obj.mask);caxis(liveArgs.cAxis);title(num2str(index));drawnow;
                    end
                    if dataYes
                        dataStorage{ss} = dat;
                    end
                    if itotYes
                        itot(ss) = sum(sum(~mask.*dat));
                    end
                    if imaxYes
                        imax(ss) = max(max(~mask.*dat));
                    end
                    if fluoYes
                        fluoData = obj.data.read(index);
                        % iterate over all counters
                        for kk = 1:length(obj.fluoCounters)
                            currCounter = obj.fluoCounters{kk}(1);
                            lower = obj.fluoCounters{kk}(2);
                            upper = obj.fluoCounters{kk}(3);
                            fluoResult.(currCounter{1})(jj) = sum(fluoData(lower{1}:upper{1},1));
                        end
                        fluoResult.sumSpectrum = fluoResult.sumSpectrum + fluoData;
                    end
                    if maxProjYes
                        maxproj = max(dat,maxproj);
                    end
                    if stxmYes
                        for ii = 1:size(stxmArgs.sel,3)
                            df(ss,ii) = stxm(dat,stxmArgs.mask,stxmArgs.sel(:,:,ii));
                            % also available as mex function
                        end
                    end
                    if dpcYes
                        [dpcx(ss),dpcy(ss)] = dpc(~mask.*dat);
                    end
                    if develYes
                        devel(ss) = NaN; %sum(sum(((~obj.mask).*dat) > 50));
                    end
                    if crystalYes
                        for ii = 1:size(crystalArgs.sel,3)
                            crystal(ss,ii) = crystals(dat-crystalArgs.bgr,crystalArgs.mask,crystalArgs.sel(:,:,ii),crystalArgs.threshold);
                        end
                    end
                    if pcaYes
                        for ii = 1:size(pcaArgs.sel,3)
                            res = pca(dat,pcaArgs.qy,pcaArgs.qz,pcaArgs.mask,pcaArgs.sel(:,:,ii));
                            w(ss,ii) = res.w;
                            phi(ss,ii) = res.angle;
                            % also available as mex function
                        end
                    end
                    if circDistYes
                        for ii = 1:size(circDistArgs.sel,3)
                            % radial averaging
%                             tmp = b1d(dat,circDistArgs.mask,circDistArgs.sel(:,:,ii),circDistArgs.grid,circDistArgs.bins);
                            % faster
                            tmp = b1d_fast(dat,circDistArgs.mask,circDistArgs.sel(:,:,ii),sortIndex_rad,n_rad,bin_rad);
                            tmp.qr = ph; % for backwards compatibility
                            
                            % remove background
                            if ~isempty(circDistArgs.bgr)
                                tmp.dat_1d = tmp.dat_1d - circDistArgs.bgr(:,ii);
                            end
                            % find correct mean and variance of circular
                            % (phase) signal
                            [circAngleTmp, circWtmp] = circular_dist(tmp.qr, tmp.dat_1d);
                            w(ss,ii) = circWtmp;
                            phi(ss,ii) = circAngleTmp;
                        end
                    end
                    if sumYes
                        sumResult = sumResult + dat;
                    end
                    if avgPerLineYes
                        avgPerLine{line} = avgPerLine{line} + dat;
                    end
                    if symmetryYes
                         [ln, symmetryBinCenter] = symmetry(dat, symmetryArgs.bgr, symmetryArgs.sigma, symmetryArgs.noise_level, obj);
                         Ln(ss,:) = ln;
                    end
                    if b1dYes
                        
                        if ~isempty(b1dArgs.bgr)
                            dat_b1d = dat - b1dArgs.bgr;
                        else
                            dat_b1d = dat;
                        end
                        
                        for ii = 1:size(b1dArgs.sel,3)
                            if size(b1dArgs.sel,3) > 1 
                                % slightly slower method
                                tmp = b1d(dat,b1dArgs.mask,b1dArgs.sel(:,:,ii),b1dArgs.grid,b1dArgs.bins);
                                b1d_res(ss,ii).dat_1d = tmp.dat_1d;
                                b1d_res(ss,ii).qr = tmp.qr;
                                b1d_res(ss,ii).x = tmp.qr;
                            else
                                % faster method
                                tmp = b1d_fast(dat_b1d,b1dArgs.mask,b1dArgs.sel(:,:,ii),sortIndex_azim,n_azim,bin_azim);
                                b1d_res(ss,ii).dat_1d = tmp.dat_1d;
                                b1d_res(ss,ii).qr = qr;
                                b1d_res(ss,ii).x = qr;                                
                            end

                        end
                    end
                    if pyfaiYes
                        tmp = pyfai(dat,pin.Results.pyFAImask,obj.pby,obj.pbz,pfsettings);
                        pyfai_res(ss).error = tmp.error;
                        pyfai_res(ss).dat_1d = tmp.dat_1d;
                        pyfai_res(ss).qr = obj.q_of_n(tmp.r);
                        if strcmp(pfsettings.to2D,'on')
                            pyfai_res(ss).dat_2d = tmp.dat_2d;
                        end
                    end
                        
                end
                fprintf(1,'\b\b\b\b\b\b\b');
            else 
                % here parallel version
                
%                 f = obj.data.compile_static_function;
%                 g = obj.process_data_static;
                
                parfor jj = 1:length(indexlist)

%                     index = indexlist(jj);
%                     fn = fnrs(index);
%                  
%                     % always read data
%                     data_temp = f(fn);
% 
%                     if strcmp(emptySub,'on')
%                         dat = mask_sliced.*(g(data_temp) - empty);
%                     else                            
%                         dat = mask_sliced.*g(data_temp);
%                     end
                    
                    

                    index = indexlist(jj);
                    fn = fnrs(index);
                    
%                     if index > (line+1)*numel(ysampling) % line done
%                         fprintf(1,'Processed line: %7d / %7d\n', line, numel(zsampling));
%                         line = line + 1;
%                     end
                    
                    % read data if necessary/requested
                    if ~contains(method,'none')

                        dat = obj.process(obj.data.read(fn));
                        % dat = obj.read_wait(fn);
                        
                        % perform pre-processing
                        if emptyYes
                            dat = dat - empty;
                        end                
                        
                        if strcmp(pin.Results.correction,'on')
                            dat = dat.*corr;
                        end
                        
                        % always apply mask
                        % dat(mask) = 0; % is this still needed?
                    end
                    
                    
                    if stxmYes
                        if size(stxmArgs.sel,3) == 1
                            df(jj) = stxm(dat,stxmArgs.mask,stxmArgs.sel(:,:,1));
                            % also available as mex function
                        end
                    end
%                     if crystalYes
%                         if size(crystalArgs.sel,3) == 1
%                             crystal(jj) = crystals(dat-crystalArgs.bgr,crystalArgs.mask,crystalArgs.sel(:,:,1),crystalArgs.threshold);
%                         end
%                     end
%                     if circDistYes
%                         if size(circDistArgs.sel,3) == 1
%                             tmp = b1d(dat,circDistArgs.mask,circDistArgs.sel(:,:,1),circDistArgs.grid,circDistArgs.bins);
%                             if ~isempty(circDistArgs.bgr)
%                                 tmp.dat_1d = tmp.dat_1d - circDistArgs.bgr(:,1);
%                             end
%                             [circAngleTmp, circWtmp] = circular_dist(tmp.qr, tmp.dat_1d);
%                             w(jj) = circWtmp;
%                             phi(jj) = circAngleTmp;
%                         end
%                     end
%                     if sumYes
%                         sumResult = sumResult + dat;
%                     end
%                     if b1dYes
%                         if size(b1dArgs.sel,3) == 1
%                             tmp = b1d(dat,b1dArgs.mask,b1dArgs.sel(:,:,1),b1dArgs.grid,b1dArgs.bins);
%                             b1d_res(jj).dat_1d = tmp.dat_1d;
%                             b1d_res(jj).qr = tmp.qr;
%                         end
%                     end
                    
                end                
            end
            
            toc;
           
            % save result
            if dataYes
                result.data = sort_scan(dataStorage,nr_rows,nr_cols,scanmode);
            end
            if maxProjYes
                result.maxProjection = maxproj; % legacy
                result.mp = maxproj;
            end
            if imaxYes
                result.imax = sort_scan(imax,nr_rows,nr_cols,scanmode);
            end
            if itotYes
                result.itot = sort_scan(itot,nr_rows,nr_cols,scanmode);
            end
            if crystalYes
                result.crystal = sort_scan(crystal,nr_rows,nr_cols,scanmode);
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
                result.pca.angle = sort_scan(phi,nr_rows,nr_cols,scanmode);
            end
            if circDistYes
                result.circDist.w = sort_scan(w,nr_rows,nr_cols,scanmode);
                result.circDist.angle = sort_scan(phi,nr_rows,nr_cols,scanmode);
            end
            if fluoYes
                fluoResult.averageSpectrum = fluoResult.sumSpectrum / length(indexlist);
                % sort every fluo counter
                for kk = 1:length(obj.fluoCounters)
                    currCounter = obj.fluoCounters{kk}(1);
                    fluoResult.(currCounter{1}) = sort_scan(fluoResult.(currCounter{1}),nr_rows,nr_cols,scanmode);
                end
                result.fluo = fluoResult;
            end
            if avgPerLineYes
                nrFramesPerLine = sum(roi,2);
                for l = 1:numel(avgPerLine)
                    avgPerLine{l} = avgPerLine{l} ./ nrFramesPerLine(l);
                end
                result.avgPerLine.bgr = avgPerLine;
                result.avgPerLine.nrFramesPerLine = nrFramesPerLine;
                result.avgPerLine.roi = roi;
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
                result.symmetry.Ln = sort_scan(Ln,nr_rows,nr_cols,scanmode);
                result.symmetry.binCenter = symmetryBinCenter;
                result.symmetry.lnHist = squeeze(sum(sum(result.symmetry.Ln,1),2));
            end
            if b1dYes
               result.b1d = sort_scan(b1d_res,nr_rows,nr_cols,scanmode);
            end
            if pyfaiYes
               result.pyfai = sort_scan(pyfai_res,nr_rows,nr_cols,scanmode);
            end
        end        
        
        
        function [stripped, varargout] = strip(obj, data, winSize, iterations, varargin)
            % STRIP  Model-indepentend background stripping for
            % one-dimensional graphs with monotone background and
            % superimposted, sharp features such as bragg reflections or 
            % spectral lines.
            %   
            %   stripped = strip(data,windowSize,iterations,qr)
            %
            % The following arguments are supported:
            %   data:: [] (required)
            %       A one-dimensional array of data points.
            %       See help: nanodiffraction.baseline.
            %
            %   windowSize:: [] (required)
            %       See help: nanodiffraction.baseline.
            %
            %   iterations:: [] (required) 
            %       See help: nanodiffraction.baseline.
            %
            %   qr:: [] (optional)
            %       Detector coordinates in terms of the radial wavevector
            %       transfer. A matrix of size [Nz,Ny] is expected.
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %       stripped::
            %           stripped background
            %
            %       stripped_2d:: (optional)
            %           2d stripped background of size [Nz,Ny]. This is
            %           only output if a two-dimensional qr-grid was given
            %           to the function as optional argument.
            %
            
            
            % pre-process
            % implement some filter here
            
            % strip
            stripped = baseline(...
                       baseline(data,winSize,iterations)...
                       ,1,500);
            
            % postprocess
            if nargin == 5
                % make 2d background
                qr = varargin{1};
                varargout{1} = obj.azimuthal_spread(stripped,qr);
            end
        end
        
        
        function [healed, il] = heal(obj, data, mask)
            % HEAL  uses make_mask_symmetric to find the correct mapping
            % from valid to invalid pixels. Takes in an n x m data matrix
            % and an n x m logical mask (1: corrupt pixel, 0: otherwise)
            % and mirrors data that is lost due to e.g. modular gaps from 
            % the point-symmetric position in the data. This function was
            % inspired by the article "Healing X-ray scattering images" by 
            % Liu et al., IUCrJ 2017, Vol. 4, Issue 4, pp. 455 - 465.
            %   
            %   healed = heal(data,mask)
            %
            % The following arguments are supported:
            %       data:: [] (required)
            %          An n x m data matrix.
            %                
            %       mask:: [] (required)
            %          An n x m logical matrix indicating bad pixels.
            %          1: invalid pixel, 0: valid pixel.
            %
            % Example:
            %   e = nanodiffraction();
            %   frame = e.read(1);
            %   mask = (frame == max(frame(:)));
            %   healed_frame = e.heal(frame,mask);
            %
            % Output arguments:
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
        
        
%         function data = read_wait(obj,fn)
%             % READ_WAIT  Useful for online data analysis when not all data
%             % is immediately ready for processing. Similar to read() but
%             % uses (nanodiffraction).wait_it and
%             % (nanodiffraction).wait_maxit and (nanodiffraction).wait to
%             % decide wether or not to wait (wait) and for how long 
%             % (wait_maxit). wait_it is the current waiting index.
%             %   
%             %   DATA = READ_WAIT(FN)
%             %
%             %   Input and output arguments as in read().
%             
%             
%             % should I wait?
%             if obj.wait && (obj.wait_it < obj.wait_maxit)
%                 try
%                     data = obj.process(obj.data.read(fn)); 
%                     obj.wait_it = 0;
%                 catch err
%                     obj.wait_it = obj.wait_it + 1;
%                     pause(1);
%                     data = obj.read(fn);
%                 end 
%             else
%                 data = obj.process(obj.data.read(fn)); 
%             end
%         end
        
        
        function [composite, parameters] = calculate_composite(obj,varargin)
            % CALCULATE_COMPOSITE  calculates a composite image from a
            % series of diffraction patterns.
            %   
            %   [composite, parameters] = calculate_composite(options) 
            %
            % The following arguments are supported:
            %     options:: [see default values below] (optional)
            %       Structure that can contain the following fields. Note,
            %       that all fields are optional. Default values that are
            %       otherwise used are as usual given in angle brackets:
            %
            %       - yCrop:: [1 obj.scan.SNy]
            %           Description.
            %
            %       - zCrop:: [1 obj.scan.SNz]
            %           Description.
            %
            %       - ySkip:: [1]
            %           Description.
            %
            %       - zSkip:: [1]
            %           Description.
            %
            % Example:
            %   e = nanodiffraction(<some settings here>);
            %   e.set_scan_info(<some settings here>);
            %   [comp,p] = e.calculate_composite(...
            %                   struct('ySkip',4,'zSkip',4));
            %
            % Output arguments:
            %       composite::
            %           Processed composite matrix.
            %
            %       parameters::
            %           Parameters used to combine the input data. Is used
            %           in combination with (display_class).composite. See
            %           help display.composite for more information
            %
            
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
            if strcmp(opts.heal,'off')
                result = obj.analyze_scan('method','data','yCrop',opts.yCrop,'zCrop',opts.zCrop,'ySkip',opts.ySkip,'zSkip',opts.zSkip,'emptySub',opts.emptySub,'correction',opts.correction);
            else
                result = obj.analyze_scan('method','data+heal','yCrop',opts.yCrop,'zCrop',opts.zCrop,'ySkip',opts.ySkip,'zSkip',opts.zSkip,'emptySub',opts.emptySub,'empty',opts.empty,'correction',opts.correction,'healmask',opts.healmask);
            end
            
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
            
%         function [bgr] = median_background(obj)
%             % MEDIAN_BACKGROUND  Calculates the median of the first lines
%             % in a scan. In some cases, this can be used as a simple
%             % estimate for a background for scattering patterns.
%             %   
%             %   [BACKGROUND] = MEDIAN_BACKGROUND() 
%             %
%             %   Output arguments:
%             %       BACKGROUND::
%             %           Calculated median background.
%             
%             
%             % line median
%             fprintf(1,'Averaging %d lines for median calculation...', ceil(50/double(obj.scan.SNy)));
%             rawdata = obj.analyze_scan('zCrop',[1 ceil(50/double(obj.scan.SNy))],'method','data');
% 
%             % reshape
%             [m,n] = size(cell2mat(rawdata.data(1)));
%             stack = zeros(m,n,numel(rawdata.data));
%             for ii = 1:numel(rawdata.data)
%                 stack(:,:,ii) = cell2mat(rawdata.data(ii));
%             end
% 
%             % median-filtered background
%             bgr = median(stack,3);
%         end
        
        
%         function robustFit = robust_fitting(obj,dat,varargin)
%             % ROBUST_FITTING  in principle, this function performs an 
%             % ordinary weighted linear least square analysis. Since it is
%             % linear, it can be written in a matrix equation form. 
%             % For more information, see:
%             % https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/06LeastSquares/general/
%             % http://de.mathworks.com/help/curvefit/least-squares-fitting.html
%             % http://www.dsplog.com/2012/02/05/weighted-least-squares-and-locally-weighted-linear-regression/
%             % 
%             %   [ROBUST_FIT] = ROBUST_FITTING(DATA,OPTS)
%             %
%             %   The following options are supported:
%             %
%             %     DATA:: []
%             %       A struct containing the fields 'dat_1d' and 'qr' is
%             %       required.
%             %
%             %     OPTS:: [] (optional)
%             %       structure that can contain the following fields:
%             %           - 'ncoeff': number of coefficients
%             %           - 'win': Window for data selection, a 2-element vector.
%             %
%             %   Output arguments:
%             %    
%             %       ROBUST_FIT::
%             %           struct containing
%             %           - f: fit amplitude
%             %           - x: fit abscissa
%             %           - c: fitted coefficients
%             
%             % parse variable arguments
%             defaults = struct('ncoeff',9,'win',[1 numel(dat.qr)],'coeff_start',0);
%             fields = fieldnames(defaults);
%             if nargin > 2
%                 opts = varargin{1};
%                 for f = 1:numel(fields)
%                     if ~isfield(opts,fields{f})
%                         opts.(fields{f}) = defaults.(fields{f});
%                     end
%                 end
%             else
%                 opts = defaults;
%             end
%             
%             n_coeff = opts.ncoeff;
%             win = opts.win;
% 
%             % data
%             x = dat.qr(win(1):win(2));
%             y = dat.dat_1d(win(1):win(2));
%             
%             % fit
%             [f,c] = lin_fit(x,y, n_coeff, @(y) 1./y.^2, opts.coeff_start);
% 
%             % save
%             robustFit.f = f;
%             robustFit.x = x;
%             robustFit.c = c;
%         end
        
        
%         function [res] = cluster(~,c,max_clusters,varargin)
%             % CLUSTER  This function starts with k random starting points 
%             % in the nx x ny x c_i - space 
%             % (e.g. a scan with 51x51 scan points and 9 fitting
%             % coefficients gives a 51x51x9 volume). It then searches for
%             % k means within this volume by minimizing the euklidian norm.
%             % 
%             %   [CLUSTERS] = CLUSTER(C,MAX_CLUSTERS,LOGSCALE)
%             %
%             %   The following options are supported:
%             %
%             %     C:: []
%             %       An n x m x k - data matrix where k = 1:K and K is the 
%             %       number of clusters. n and m are row and col position in
%             %       the scan.
%             %
%             %     MAX_CLUSTERS:: [] 
%             %       Maximum numbers of clusters (i.e. means).
%             %
%             %     LOGSCALE:: [1] (optional)
%             %       If set to 1, k-means will be displayed in log-scale.
%             %
%             %   Output arguments:
%             %    
%             %       See help cluster_analysis.
%             
%             if nargin == 4
%                 logscale = varargin{1};
%             else
%                 logscale = 1;
%             end
%             
%             nx = size(c,2);
%             ny = size(c,1);
%             nz = size(c,3);
%             
%             % bring c into form of (coeff,datapoint)
%             c_reshaped = zeros(nz,nx*ny);
%             for ii = 1:size(c,3)
%                 c_reshaped(ii,:) = squeeze(reshape(c(:,:,ii),[1,nx*ny]));
%             end
%             
%             res = cluster_analysis(c_reshaped,max_clusters,logscale);
%             res.clusters = reshape(res.cluster,ny,nx);
%             
%         end
        
        function mask = process_mask(obj,mask)
            % PROCESS_MASK  Used for processing logical masks. The mask
            % will be adapted to the current setting of a detector ROI and
            % binning factors.
            %
            %   processed_mask = process_mask(mask)
            %
            % The following arguments are supported:
            %   mask:: [] (required)
            %       A two-dimensional mask. The dimensions of the array
            %       should correspond to the detector size.
            %   
            % Example:
            %   e = nanodiffraction();
            %   mask = e.read(1) > 1000;
            %   e.set_roi_and_binning('binning','on','biny',4,'binz',4);
            %   new_mask = e.process_mask(mask);
            %
            % Output arguments:
            %   processed_mask:: The cropped and rebinned logical mask.
            %   
            
            
            mask = obj.process(mask);
            mask = min(mask,1);
        end
            
        
        function sel_out = set_selection(obj, varargin)
            % SET_SELECTION  The function accepts as many selections as 
            % needed as arguments, as long as their dimensions are 
            % identical. 
            % Note, that a warning will be shown, if the size of the mask 
            % does not correspond to the current detector dimensions. 
            % 
            %   [sel_out] = set_selection(sel, sel, sel, ...)
            %
            % The following arguments are supported:
            %     sel:: [] (required)
            %       A two-dimensional array. The input should be a logical 
            %       array. 
            %
            %     sel:: [] (optional)
            %       Additional selections can be given. All selections will 
            %       then be logically combined (OR combination [|]).
            %
            %     sel:: [] (optional)
            %       The amount of selections is not limited, all selections 
            %       will be combined to one.
            %
            % Example:
            %   e = nanodiffraction(<some settings>);
            %   sel_small = e.radial_mask('r1',0,'r2',1.2,'grid',e.qr);
            %   sel_large = e.radial_mask('r1',2.4,'r2',3.2,'grid',e.qr);
            %   e.set_selection(sel_small,sel_large);
            %
            % Output arguments:
            %       sel_out::
            %           Current selection.
            %
            % Note, that selections are stored in:
            %   (nanodiffraction_class).sels
            %
            % Note, that more selections can be used, for example, when
            % performing a pca analysis or stxm. Additional selections
            % can be added by calling 
            %	(nanodiffraction).add_selection(sel)
            %
            % For more help, please read:
            %   help nanodiffraction.add_selection
            
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
        
        
        function selections = add_selection(obj, sel)
            % ADD_SELECTION  The function adds a logical selection matrix 
            % to the current stack of selections.
            % 
            %   selections = add_selection(selection)
            %
            % The following arguments are supported:
            %     selection:: []
            %       A two-dimensional array. The input should be a logical
            %       array. 
            %
            % Example:
            %   e = nanodiffraction(<some settings>);
            %   sel_small = e.radial_mask('r1',0,'r2',1.2,'grid',e.qr);
            %   sel_large = e.radial_mask('r1',2.4,'r2',3.2,'grid',e.qr);
            %   sel_all = ones(e.Nz,e.Ny);
            %   e.set_selection(sel_all);
            %   e.add_selection(sel_small);
            %   e.add_selection(sel_large);
            %
            % Output arguments:
            %       selections::
            %           Current set of selections.
            %
            %   Note, that the selections are stored in:
            %       (nanodiffraction_class).sels
            
            % add selection if dimensions fit to previously added selection
            if numel(sel) == numel(obj.sels(:,:,1))
                obj.sels = cat(3,obj.sels,sel);
            else
                error('Selection matrix does not have the same size as the selections stored in obj.sels');
            end

            selections = obj.sels;
            
            if max(size(sel) == [obj.Nz_orig,obj.Ny_orig])
                obj.sels_orig = obj.sels;     
            end
            
        end
        
                    
        function mask_out = set_mask(obj, varargin)
            % SET_MASK  The function accepts as many masks as needed as
            % arguments, as long as their dimensions are identical. 
            % Note, that a warning will be shown, if the size of the mask 
            % does not correspond to the current detector dimensions. 
            % 
            %   [mask_out] = set_mask(mask, mask, mask, ...)
            %
            % The following arguments are supported:
            %     mask:: [] (required)
            %       A two-dimensional array.The input should be a logical
            %       array. 
            %
            %     mask:: [] (optional)
            %       Additional masks can be given. All masks will then be 
            %       logically combined (OR combination [|]).
            %
            %     mask:: [] (optional)
            %       The amount of masks is not limited, all masks will then
            %       be logically combined (OR combination [|]).
            %
            % Example:
            %   e = nanodiffraction(<some settings>);
            %   mask1 = e.radial_mask('r1',0,'r2',0.2,'grid',e.qr);
            %   mask2 =
            %   e.azimuthal_mask('phi1',30,'phi2',56,'grid',e.phi2);
            %   e.set_mask(mask1,mask2);
            %
            % Output arguments:
            %       mask_out::
            %           Combined mask.
            %
            % Note, that the detector mask is stored in:
            %       (nanodiffraction_class).mask
            %
            
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
        
        
        function corr_out = set_corr(obj, corr)
            % SET_CORR  Defines a correction matrix to be applied in a
            % subsequent analysis.
            % Note, that a warning will be shown, if the size of the
            % correction matrix does not correspond to the current detector 
            % dimensions. 
            % 
            %   [corr_out] = set_corr(corr)
            %
            % The following arguments are supported:
            %     corr:: [] 
            %       A two-dimensional array. Input can a double array.
            %
            % Example:
            %   e = nanodiffraction(<some settings>);
            %   e.set_scan_info('SNy',512,'SNz',651);
            %   mask = e.masks | e.azimuthal_mask('phi1',-20,'phi2',20);
            %   bs = e.analyze_scan('method','average','zCrop',[1 2],'yCrop',[2 1079]);
            %   rad_bgr = b1d(bs.avg,mask,[],e.qr,512);
            %   vq2d = e.azimuthal_spread(rad_bgr.dat_1d,rad_bgr.qr);
            %   e.set_corr(~e.masks.*(vq2d./bs.avg));
            %
            % Output arguments:
            %       corr_out::
            %           Active correction matrix
            %
            % Note, that the correction matrix is stored in: 
            %   (nanodiffraction_class).corr
            % after set_corr was used.
            
 
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
            % SET_ROI_AND_BINNING  All following functions require a standard 
            % pixel mask that mask hot pixels and intermodular gaps of 
            % detectors. In addition, a calibration could be added that 
            % corrects for solid angle deviations at larger scattering 
            % angles.
            %
            %   set_roi_and_binning(pair_values_argument_list)
            %
            % The following arguments are accepted:
            %      detCalibration:: [1] (optional)
            %        mask taking into account corrections e.g. for solid 
            %        angles
            %
            %      detectorRoi:: [off] (optional)
            %        Can be toggled between on and off
            %
            %      detRoiY:: [start end] (optional)
            %        [start end] in pixel units
            %
            %      detRoiZ:: [start end] (optional)
            %        [start end] in pixel units  
            %
            %      binning:: [off] (optional)
            %        Can be toggled between on and off
            %
            %      binx:: [1] (optional)
            %        Binning ratio along y
            %
            %      biny:: [1] (optional)
            %        Binning ratio along z
            %
            % Example:
            %   e = nanodiffraction();
            %   e.set_roi_and_binning('detectorRoi','on',...
            %       'detRoiY',around_y(100),...
            %       'detRoiZ',around_z(100),...
            %       'binning','on',...
            %       'biny',2,...
            %       'binz',2);
            %
            % Output arguments:
            %   This function does not return any arguments.
            %
        
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
            disp('Resetting number of pixels, pixel size and primary beam coordinates to original values');
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
            % analyses. For p10, one could use spec_get_scan_info, 
            % therefore choose p10 for the beamline argument. 
            % If one chooses to do so, please also provide a "scanNo",
            % "newfile" and "specpath" argument, such that 
            % "spec_get_scan_info" can be used in turn.
            % However, this is usually too complicated and the beamline
            % argument can be disregarded completely. 
            %
            %   set_scan_info(pair_values_argument_list)
            %
            % The following arguments are supported:
            %       SNy, SNz:: [1] (required)
            %          number of scan points along fast and slow axis
            %   
            %       stepy, stepz:: [1e-6] (optional)
            %          Stepsize along both axes in meters
            %
            %       firstFile:: [1] (optional)
            %          In some special cases, the first file of a scan is
            %          not 1. If this is the case, please change the
            %          "firstFile" argument accordingly (however, this 
            %          should not normally happen)
            %
            % Example:
            %   e = nanodiffraction();
            %   e.set_scan_info('SNy',51,'SNz',51);
            %
            % Output arguments:
            %   This function does not return any arguments.
            %
            
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
            %   empty = calculate_empty(pair_values_argument_list)
            %
            % The following arguments are supported:
            %     roi:: [] (optional)
            %        Logical mask that defines a region-of-interest
            %
            %      yCrop:: [1 NR_OF_SCAN_POINTS_ALONG_FAST_AXIS] (optional)
            %        [start end] in pixel units
            %
            %      zCrop:: [1 NR_OF_SCAN_POINTS_ALONG_SLOW_AXIS] (optional)
            %        [start end] in pixel units  
            %
            %      ySkip:: [1] (optional)
            %        only every n-th (n integer) scan point along fast axis
            %        shall be evaluated
            %
            %      zSkip:: [1] (optional)
            %        only every n-th (n integer) scan point along slow axis
            %        shall be evaluated
            %
            % Example:
            %   e = nanodiffraction();
            %   empty = e.calculate_empty('yCrop',[1 10],'zCrop',[1 10]);
            %
            % Output arguments:
            %   This function does not return any arguments.
            %
        
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
            obj.empty = empty;
        end
        
       
        function dat = process(obj,dat)
            % PROCESS  Used for processing detector images. The frame
            % will be adapted to the current setting of a detector ROI and
            % binning factors.
            %
            %   processed_data = process(data)
            %
            % The following arguments are supported:
            %   data:: [] (required)
            %       A two-dimensional detector image. The dimensions of the 
            %       array should correspond to the detector size.
            %   
            % Example:
            %   e = nanodiffraction();
            %   frame = e.read(1);
            %   e.set_roi_and_binning('binning','on','biny',4,'binz',4);
            %   processed_frame = e.process(frame);
            %
            % Output arguments:
            %   processed_data:: The cropped and rebinned detector image.
            %   
            
            
            if obj.detRoiZ(1) > size(dat,1) || obj.detRoiZ(2) > size(dat,1) || ...
                    obj.detRoiY(1) > size(dat,2) || obj.detRoiY(2) > size(dat,2)
                fprintf('Size of matrix: %d (row) x %d (col)\n', size(dat,1), size(dat,2));
                fprintf('Roi (y): [%d %d], Roi (z): [%d %d]\n', obj.detRoiY(1),obj.detRoiY(2),obj.detRoiZ(1),obj.detRoiZ(2));
                error('Sizes do not match');
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
        
        
%         function g = process_data_static(obj)
%             % PROCESS_DATA_STATIC  slices all variables for use in parfor
%             % loop (under development)
%             %
%             %   FUNCTIONHANDLE = PROCESS_DATA_STATIC()
%             
%             biny = obj.biny;
%             binx = obj.binx;
%             detRoiZ = obj.detRoiZ;
%             detRoiY = obj.detRoiY;
%             if strcmp(obj.detectorRoi,'on') && ~strcmp(obj.binning,'on') 
%                 g = @(dat) dat(detRoiZ(1):detRoiZ(2),detRoiY(1):detRoiY(2));
%                 return;
%             end
%             if strcmp(obj.binning,'on') && ~strcmp(obj.detectorRoi,'on')
%                 [m,n] = deal(obj.Nz_orig,obj.Ny_orig);
%                 nnew = n - mod(n,obj.binx);
%                 mnew = m - mod(m,obj.biny);
%                 g = @(dat) transpose(reshape(sum(reshape(transpose(reshape(sum( reshape(dat(1:mnew,1:nnew),biny,[]),1),mnew/biny,[])),binx,[]),1),nnew/binx,[]));
%                 return;
%             end
%             if strcmp(obj.binning,'on') && strcmp(obj.detectorRoi,'on')
% %                 [m,n] = deal(obj.Nz,obj.Ny);
%                 [m,n] = deal(detRoiZ(2)-detRoiZ(1)+1,detRoiY(2)-detRoiY(1)+1);
%                 
%                 nn = n - mod(n,obj.binx);
%                 mm = m - mod(m,obj.biny);
% 
%                 % now it is divisible by binx, biny
%                 g = @(dat) reshape(sum(sum(...
%                     reshape(...
%                         dat(detRoiZ(1):(detRoiZ(2) - mod(n,obj.binx)),detRoiY(1):(detRoiY(2) - mod(n,obj.binx))),...
%                         [biny mm/biny binx nn/binx]),...
%                     1),3),[mm/biny nn/binx]);
% 
%                 return;
%             else
%                 g = @(dat) dat;
%                 return;
%             end
%         end
        
        
        function mask = radial_mask(obj,varargin)
            % RADIAL_MASK  defines a radial mask (e.g. useful for a PCA
            % analysis). Creates a circular or ring-shaped mask and stores 
            % it in (nanodiffraction_class).sels.
            %
            %   mask = radial_mask(pair_values_argument_list) 
            %
            % The following arguments are accepted:
            %
            %       r1:: [-1] (optional)
            %          Inner radius
            %
            %       r2:: [232] (optional)
            %          Outer radius (Inverse radius)
            %
            %       grid:: [nanodiffraction_class.qr) (optional)
            %           Two-dimensional detector coordinates.
            %
            %       A short note regarding the usage of 'r1' and 'r2': 
            %       If both parameters are given, the function will 
            %       generate logical 1, if r is within r1 and r2
            %       If only r1 is given, logical mask will be 1 for r < r1
            %       If only r2 is given, logical mask will be 1 for r > r2
            %
            % Example:
            %   e = nanodiffraction(<some settings>);
            %   mask = e.radial_mask('r1',0.2,'r2',1.2,'grid',e.qr);
            %
            % Output arguments:
            %   mask:: A two-dimensional logical radial mask.
            %   
        
            % parse input
            pin = inputParser;
            addOptional(pin,'grid',obj.R);
            addOptional(pin,'r1',-1);
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
            % AZIMUTHAL_MASK  defines an azimuthal mask (e.g. useful for 
            %   PCA analysis). This function creates a circular or 
            %   ring-shaped mask. Note, that the angles in degrees follow 
            %   the mathematical convention (counter-clockwise and 0° at
            %   3 o'clock).
            %
            %   mask = azimuthal_mask(pair_values_argument_list) 
            %
            %   The following parameters are accepted:
            %
            %       phi1:: [65] (optional)
            %          Start angle in degrees.
            %
            %       phi2:: [232] (optional)
            %          End angle in degrees
            %
            %       grid:: [nanodiffraction_class.phi2) (optional)
            %           Two-dimensional detector coordinates.            
            %
            %       A short note regarding the usage of 'r1' and 'r2': 
            %       If both parameters are given, the function will 
            %       generate logical 1, if phi is within phi1 and phi2
            %       If only phi1 is given, the logical mask will be 1 for 
            %       phi < phi1
            %       If only phi2 is given, the logical mask will be 1 for 
            %       phi > phi2
            %
            % Example:
            %   e = nanodiffraction(<some settings>);
            %   mask = e.azimuthal_mask('phi1',30,'phi2',45,'grid',e.phi2);
            %
            % Output arguments:
            %   mask:: A two-dimensional logical azimuthal mask.
            %               
        
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
        
                
%         function [cm] = estimate_primary_beam(obj,data)
%             % ESTIMATE_PRIMARY_BEAM  function that estimates primary beam 
%             % position based on center of mass (under development)
%         
%             % estimate primary beam position
%             [obj.pbz,obj.pby] = CM_nano(data); % output format [col, row]
% 
%             % display result if needed
%             prompt = 'Do you wish to see the result? (y/n)';
%             answer = input(prompt,'s');
%             if strcmp(answer,'y')
%                 imagesc(log10(abs(data)));
%             end
%         end
        
        
%         function [orientation] = orientation(~,pca)
%             % ORIENTATION  shorthand to get the orientation angle instead 
%             % of the x- and y-component of the angle.
%             %
%             %   PHI = ORIENTATION(PCA) accepts a struct that is formatted
%             %   as the output of analyze scan, as 
%             %       pca - output1
%             %           - output2...
%             
%             orientation = atan(pca.ycomp./pca.xcomp);
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
            %   linear_index = sub2ind(col,row) 
            %
            % The following arguments are accepted:
            %       row:: [] (required)
            %          row index in scan (z-axis).
            %
            %       col:: [] (required)
            %          column index in scan (y-axis).
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   linear_index:: Linear index.
            %
            
            roi = zeros(obj.scan.SNz,obj.scan.SNy);
            roi(z,y) = 1;
            indexlist = find(flipud(rot90(roi)));
            ind = indexlist(1);
        end
        
        
        function [row,col] = ind2sub(obj,ind,sny)
            % IND2SUB  calculates the row and column index based on the
            % linear index in the current scan.
            %
            %   [row,col] = ind2sub(linear_index, sny)
            %
            % The following arguments are supported:
            %
            %       linear_index:: [] (required)
            %          Linear index. The first index in the scan is 1.
            %
            %       sny:: [] (required)
            %          Number of scan points per line (fast axis).
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   row:: Row index in scan.
            %   
            %   col:: Column index in scan.
            %
            row = floor(double(ind-1)/double(sny))+1; % 1 -> SNz
            col = mod((ind-1),double(sny))+1; % 1 -> SNy
        end
        
        
        function [parameter] = moment_from_selection(obj,varargin)
            % MOMENT_FROM_SELECTION  In an active figure, a roi can be
            % selected and the mean, standard deviation, variance and the
            % median can be obtained.
            %
            %   parameter = moment_from_selection(opts)
            %
            % The following arguments are supported:
            %
            %       opts:: ['mean']  (optional)
            %          Moment that should be calculated. By default, the
            %          mean will be calculated if opts is not set.
            %          The following options are available:
            %               'mean','std','var','median'
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   parameters:: Description missing.
            
            if nargin > 1
                method = varargin{1};
            else 
                method = 'mean';
            end
            
            % get current figure
            figure(gcf);
            axhandle = get(gca,'Children');
            data = axhandle.CData;
            
            % select roi
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
        
        
        function clicktool(obj,vis,map,varargin)
            % CLICKTOOL  Once a parameter map of the data is shown in the
            % active figure, clicktool can be called and a diffraction
            % pattern or 1d saxs curves are loaded, based on the scan point
            % that was clicked on.
            %
            %   clicktool(vis,map,optional)
            %
            % The following arguments are supported:
            %
            %       vis:: [] (required)
            %          visualization tool (required).
            %
            %       map:: [] (required)
            %          Depending on the crop and skip parameters used
            %          during the last scan analysis, the click location
            %          might have to be mapped onto the corresponding
            %          position in the scan. It map is empty, 1x1 mapping
            %          will be used, otherwise, map has to be an array
            %          containing the following parameters 
            %           [yCrop zCrop ySkip zSkip]
            %
            %       optional:: [] (optional)
            %          If no optional parameters are given, then a
            %          diffraction pattern will loaded, otherwise, an
            %          arbitrary number of structs can be given, containing
            %          an m:SNy:SNz data block (m: length of single saxs
            %          curve). The struct has to contain a data block and a
            %          scale, hence, should be passed e.g. as 
            %           struct('data',data,'scale',sc)
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   This function does not return any arguments.
            %
           
            
            % get figure of parameter map and figure of diffraction pattern
            fig_map = figure(1);
            new_fig = figure(2);
            
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
            
            
            
            % now activate parameter map figure
            figure(fig_map); 
            
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
                    
                % output the result in output figure
                figure(new_fig); drawnow;
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
                                frame_in_scan = (r-1)*obj.scan.SNy + c;
                                image_nr = obj.scan.fnrs(frame_in_scan);
                                disp(['Reading frame #' num2str(frame_in_scan)]);
                                vis.diffraction(obj.readm(image_nr),struct('process','off'));
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
                
                % switch back to old figure, note that drawnow is essential
                % here: https://de.mathworks.com/matlabcentral/answers/216085-ginput-not-working-with-r2015a-for-handling-two-figures
                figure(fig_map); drawnow;
                
            end
        end
        
        
        function data = read(obj,fn)
            % READ  reads data from the data module. Note, that data is not
            % treated in any way (binning, cropping, etc.).
            %
            %   raw_data = read(fn) reads single frame.
            %
            % The following arguments are supported:
            %       fn:: [] (required)
            %          Frame number (>0). Note that the frame number always
            %          identifies the frame within a scan. This should not
            %          be confused with the file number used in the
            %          files-module. There, the file number is associated
            %          with a specific file within a newfile-session. For
            %          this reason, the file number could start with 0.
            %
            % Example:
            %   e = nanodiffraction();
            %   first_frame = e.read(1);
            %
            % Output arguments:
            %   raw_data:: Raw detector frame.
            %
            
            data = obj.data.read(fn);
        end
        
        function data = readp(obj,fn)
            % READP  reads data from data module. Note, that data is
            % processed using the built-in function 
            % (nanodiffraction_class).process.
            %
            %   processed_data = readp(fn)
            %
            % The following arguments are supported:
            %
            %       fn:: [] (required)
            %          Frame number (>0). Note that the frame number always
            %          identifies the frame within a scan. This should not
            %          be confused with the file number used in the
            %          files-module. There, the file number is associated
            %          with a specific file within a newfile-session. For
            %          this reason, the file number could start with 0.
            %
            % Example:
            %   e = nanodiffraction();
            %   e.set_roi_and_binning('binning','on','biny',2,'binz',2);
            %   first_frame = e.readp(1);
            %
            % Output arguments:
            %   processed_data:: Processed detector frame.
            %
            
            data = obj.process(obj.data.read(fn));
        end
        
        function data = readm(obj,fn)
            % READM  reads data from data module. Note, that data is
            % processed using the built-in function 
            % (nanodiffraction_class).process and the detector mask is 
            % applied on the data frame.
            %
            %   processed_data = readm(fn) 
            %
            % The following arguments are supported:
            %
            %       fn:: [] (required)
            %          Frame number (>0). Note that the frame number always
            %          identifies the frame within a scan. This should not
            %          be confused with the file number used in the
            %          files-module. There, the file number is associated
            %          with a specific file within a newfile-session. For
            %          this reason, the file number could start with 0.
            %
            % Example:
            %   e = nanodiffraction();
            %   e.set_roi_and_binning('binning','on','biny',2,'binz',2);
            %   e.set_mask(e.read(1) > 1000);
            %   first_frame = e.readm(1);
            %
            % Output arguments:
            %   processed_data:: Processed detector frame.
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
            %   [symm_mask, mask_orig, mask, good_zero, bad_zero] = make_mask_symmetric(mask)
            %
            % The following arguments are supported:
            %       mask:: []  (required)
            %          A valid logical detector mask.
            %            
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   symm_mask:: Description missing.
            %
            %   mask_orig:: Description missing.
            %
            %   mask:: Description missing.
            %
            %   good_zero:: Description missing.
            %
            %   bad_zero:: Description missing.
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
            %   image = azimuthal_spread(data_1d,qr) 
            %
            % The following arguments are supported:
            %       dat_1d:: [] (required)
            %          One-dimensional structure factor as obtained e.g.
            %          from any azimuthal integration routine such as b1d 
            %          or pyfai.
            %
            %       qr:: [] (required)
            %           Corresponding x-axis, given in units of the radial
            %           wavevector transfer qr (reciprocal nanometers).
            %            
            % Example:
            %   e = nanodiffraction(<some settings>);
            %   e.set_scan_info('SNy',512,'SNz',651);
            %   mask    = e.masks | e.azimuthal_mask('phi1',-20,'phi2',20);
            %   bs      = e.analyze_scan('method','average','zCrop',[1 2],'yCrop',[2 1079]);
            %   rad_bgr = b1d(bs.avg,mask,[],e.qr,512);
            %   image   = e.azimuthal_spread(rad_bgr.dat_1d,rad_bgr.qr);
            %
            % Output arguments:
            %   image:: Filtered image.
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
            %   stitched = stitch(dataset,stitchy,stitchz,opts) 
            %
            % The following arguments are supported:
            %
            %       dataset:: [] (required)
            %          One-dimensional cell array data{:} that contains
            %          two-dimensional diffraction patterns of identical
            %          dimensions.
            %
            %       stitchy:: [] (required)
            %           Number of images to be concatenated horizontally.
            %
            %       stitchz:: [] (required)
            %           Number of images to be concatenated vertically.
            %            
            %       opts:: [see default values below] (optional)
            %           Structure that can contain the following fields.
            %           Note, that all fields are optional. Default values
            %           that are otherwise used are as usual given in angle
            %           brackets:
            %               
            %           - win:: []
            %           4-element vector containing the lower and upper
            %           limits (in pixel units) along the horizontal and
            %           lower and upper limits along the vertical dimension
            %           of the scattering pattern. 
            %           Example: [100 900 100 900]
            %
            %           - append:: [0]
            %           Appends N fields of zeros.
            %
            %           - prepend:: [0]
            %           Prepends N fields of zeros.
            %
            %           - replaceEmpty:: [0]
            %           If a cell array at a random location is empty, it
            %           can be replaced by a field of zeros.
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   stitched:: A stitched image.
            
            defaults = struct('win',[],'append',0,'prepend',0,'replaceEmpty',0);
            fields = fieldnames(defaults);
            if nargin > 4
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
            if opts.replaceEmpty > 0
                ii = 1;
                while isempty(data{ii})
                    ii = ii + 1;
                end
                [snz,sny,snx] = size(data{ii});
                
                for ii = 1:numel(data)
                    if isempty(data{ii})
                        data{ii} = zeros(snz,sny,snx);
                    end
                end
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
        
        
        function [im_out,corr] = remove_semitransparent_object(obj,im_in,phi1,phi2)
            % REMOVE_SEMITRANSPARENT_OBJECT  Description.
            %
            %   [image_out,correction] = remove_semitransparent_object(image,phi1,phi2) 
            %
            % The following arguments are supported:
            %
            %       image:: [] (required)
            %          Description.
            %
            %       phi1:: [] (required)
            %           Description.
            %
            %       phi2:: [] (required)
            %           Description.
            %            
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   image_out:: A corrected image.
            %
            %   correction:: A correction matrix.
            %
            
            
            % perform azimuthal average on phi-range
            sel = obj.azimuthal_mask('phi1',phi1,'phi2',phi2);
            azim_average = b1d(im_in,obj.mask,sel,obj.qr,500);
            
            % spread out 1d data on 2d grid
            spreaded = obj.azimuthal_spread(azim_average.dat_1d,azim_average.qr);
            
            % calculate correction
            corr = spreaded./im_in;
            
            % final image
            im_out = im_in .* corr;
            
        end
        
        
        function [img] = simulate_agbh(obj,varargin)
            % SIMULATE_AGBH  Description.
            %
            %   simulation = simulate_agbh(opts) 
            %
            % The following arguments are supported:
            %
            %       opts:: [see default values below] (optional)
            %          Description.
            %            
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   simulation:: A simulated diffraction pattern.
            %
            
            % d spacing of silver behenate
            d = 5.8380e-9; % meters
            sigma = 0.1e-9; % meters
            q = 2*pi/d*1e-9; % rec. nm
            dq = 0.1; % rec. nm
            
            % photon flux
            n_photons_max = 100;
            
            if nargin == 4
                one_or_two_d = '2d';
                alpha = varargin{1};
                beta = varargin{2};
                gamma = varargin{3};
            else
                one_or_two_d = '1d';
            end
            
            switch one_or_two_d
                case '1d'

                    qr_max = max(obj.qr(:));
                    n_peaks = ceil(qr_max/q);
                    

                    % no background
                    qr = linspace(0,qr_max,500);
                    I = zeros(1,numel(qr));

                    for ii = 1:n_peaks
                        q0 = ii*q;
                        Igauss = n_photons_max*exp(-((qr - q0)/dq).^2);
                        I = I + Igauss;
                    end
                    
                    img = obj.azimuthal_spread(I,qr);
                    
                case '2d'
                    qr = obj.qr_tilt(alpha,beta,gamma);
                    qr_max = max(qr(:));
                    n_peaks = ceil(qr_max/q);
                    I = zeros(obj.Nz, obj.Ny);
                    for ii = 1:n_peaks
                        q0 = ii*q;
                        Igauss = n_photons_max.*exp(-((qr - q0)/dq).^2);
                        I = I + Igauss;
                    end
                    
                    img = I;
                otherwise
                    error('nothing done');
            end
                
        end
        
        
        
        function [I] = simulate_actomyosin(obj,phi0,varargin)
            % SIMULATE_ACTOMYOSIN  Description.
            %
            %   simulation = simulate_actomyosin(opts) 
            %
            % The following arguments are supported:
            %
            %       phi0:: [see default values below] (optional)
            %          Description.
            %            
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   simulation:: A simulated diffraction pattern.
            %
            
            % d spacing of actomyosin
            d = 40e-9; % meters
            sigma = 0.034; % rec. nm
            q0 = 2*pi/d*1e-9; % rec. nm
            
            % anisotropy
            dphi = 20;
            
            % photon flux
            n_photons_max = 100;
            
            IGauss = n_photons_max.*exp(-((obj.qr - q0)/sigma).^2);
            AzimGauss = 1.*exp(-((obj.phi - phi0)/dphi).^2);
            I = IGauss.*AzimGauss;
            
            if nargin ~= 3
                I = I + obj.simulate_actomyosin(phi0+180,1);
                I = I + obj.simulate_actomyosin(phi0+360,1);
            end
                
        end
        
        
        
        function [streak] = simulate_streak(obj, varargin)
            % SIMULATE_STREAK  simulates a gaussian streak at the location
            % of the primary beam.
            %
            %   [streak] = simulate_streak(opts)
            %
            % The following arguments are supported:
            %       opts:: [see default values below] (optional)
            %           Structure that can contain the following fields.
            %           Note, that all fields are optional. Default values
            %           that are otherwise used are as usual given in angle
            %           brackets.
            %
            %           - long:: [0.005]
            %               long axis length of streak
            %
            %           - short:: [0.05] 
            %               short axis length of streak
            %
            %           - y:: []
            %               Description.
            %
            %           - z:: []
            %               Description.
            %
            %           - Ny:: []
            %               Number of horizontal detector pixels 
            %
            %           - Nz:: []
            %               Number of vertical detector pixels 
            %
            %           - angle:: [30] 
            %               rotation angle of the streak
            %
            %           - I0:: [1]
            %               number of incident photons
            %
            % Example:
            %   e = nanodiffraction();
            %   streak = e.simulate_streak();
            %   imagesc(streak); axis image;
            %
            % Output arguments:
            %   streak:: A diffraction pattern containing a cigar-shaped
            %   streak.
            %
            
            defaults = struct('long',0.005,'short',0.05,'y',obj.yaxis,'z',obj.zaxis,'Ny',obj.Ny,'Nz',obj.Nz,'angle',30,'I0',1);
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
            
            ffgauss   = @(x,a) (exp((-1/2)*(x*a).^2)); % Gauss /(sqrt(2*pi)/a)
            rotmat    = @(alpha) [[cosd(alpha) -sind(alpha)];[sind(alpha) cosd(alpha)]];
            [long,short,y,z,ny,nz,angle,I0] = split_struct(opts,fields);
            
            % generate mesh
            Y = repmat(y,nz,1);
            Z = repmat(z,ny,1)';
            
            % list of points
            v = [obj.helper.to1dcol(Y),obj.helper.to1dcol(Z)];
            v = rotmat(angle)*v';
            Y = obj.helper.to2d(v(1,:),ny,nz);
            Z = obj.helper.to2d(v(2,:),ny,nz);
            
%             % rotate mesh
%             for ii = 1:ny
%                 for jj = 1:nz
%                     tmp = rotmat(30)*transpose([Y(jj,ii),Z(jj,ii)]);
%                     Y(jj,ii) = tmp(1); Z(jj,ii) = tmp(2);
%                 end
%             end
            streak = I0 * ffgauss(Y,long).*ffgauss(Z,short);
%             imagesc(streak)
        end
        
        
        function [bgr] = simulate_airscattering(obj, varargin)
            % SIMULATE_AIRSCATTERING  simulates the effect of air
            % scattering. 
            %
            %   [bgr] = simulate_airscattering(opts)
            %
            % The following arguments are supported:
            %       opts:: [see default values below] (optional)
            %           Structure that can contain the following fields.
            %           Note, that all fields are optional. Default values
            %           that are otherwise used are as usual given in angle
            %           brackets. All distances are given in millimeters.
            %
            %           - l_pre:: [10]
            %               Pre-distance.
            %
            %           - l_post:: [10]
            %               Post-distance.
            %
            %           - l_sd:: [500]
            %               Distance sample-to-detector
            %
            %           - l_sampl:: [0.1]
            %               Sampling interval.
            %   
            %           - ap_size:: [0.07]
            %               Aperture size.
            %
            %           - I0:: [1]
            %               number of incident photons
            %
            % Example:
            %   e = nanodiffraction();
            %   airscat = e.simulate_airscatting();
            %
            % Output arguments:
            %   bgr:: Simulated air-scattering background.
            %
            
            defaults = struct(...            
                'l_pre', 10, 'l_post', 10, 'l_sd', 500, 'l_sampl', 0.1,...
                'I0', 1, 'ap_size', 0.07);
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
            
            % shorthands
            [l_pre,l_post,l_sd,l_sampl,I0,ap_size] = split_struct(opts,fields);
            
            n_slices_pre = round(l_pre / l_sampl);
            n_slices_post = round(l_post / l_sampl);
            bgr_slice = I0 / (n_slices_pre + n_slices_post);
            bgr_post = n_slices_post * bgr_slice;

            % detector
%             y_det = (-(ny-1)/2:(ny-1)/2) * px;
%             z_det = (-(nz-1)/2:(nz-1)/2) * px;
%             [y_det, z_det] = meshgrid(y_det,z_det);
%             r_det = sqrt(y_det.^2 + z_det.^2);
            r_det = obj.R;

            bgr = zeros(size(r_det));
            for i = 1:n_slices_pre
                r_cone = ap_size/2 * (l_pre - i*l_sampl + l_sd ) / (l_pre - i*l_sampl);
                bgr = bgr + (r_det < r_cone)*bgr_slice;
            end
            bgr = bgr + bgr_post;
%             bgr_saxs = b1d(bgr,[],[],r_det);

%             d = display;
%             subplot(1,2,1);d.stxm(bgr,struct('sampl',10));
%             subplot(1,2,2);d.saxs(bgr_saxs,'mode','loglog');grid on; legend('simulated air scattering');
%             set(gcf,'Position',[0 0 1400 400]);
        end
        
        
        function bs_size_det_q = simulate_bs(obj, bs_size, bs_dist)
            % SIMULATE_BS  Estimates the size of the beamstop that can be
            % expected projected on the detector.
            %
            %   [bs_size_det_q] = simulate_bs(bs_size, bs_dist)
            %
            % The following arguments are supported:
            %
            %       bs_size:: []
            %           Diameter of beamstop in (mm)
            %
            %       bs_dist:: []
            %           Distance in (mm) behind focal spot
            %
            % Example:
            %   e = nanodiffraction();
            %   bs = e.simulate_bs(0.5,100);
            %
            % Output arguments:
            %   bs_size_det_q:: Size of a circular beamstop projected onto
            %   the detector.
            
            % convert to m
            bs_size = bs_size*0.001;
            bs_dist = bs_dist*0.001;
            
            % bs size in (mm) on detector
            bs_size_det = bs_size * obj.detDistance / bs_dist;
            
            % according q
            bs_size_det_q = obj.q_of_n(bs_size_det*0.5 / obj.pixelsize);
            
            
        end
        
        function arr_out = verify_array(obj, arr_in, ref, process_handle)
            % VERIFY_ARRAY  Description.
            %
            %   output_array = verify_array(input_array, ref, process_handle)
            %
            % The following arguments are supported:
            %
            %       input_array:: []
            %           Description.
            %
            %       ref:: []
            %           Description.
            %
            %       process_handle:: []
            %           Description.
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   output_array:: Description.
            %
            
            % determine input size
            size_full = size(arr_in);
            size_slice = size_full(1:2);
            nr_slices = size(arr_in,3);
            arr_out = zeros(obj.Nz, obj.Ny,nr_slices);
            
            if min(size_slice ~= [obj.Nz,obj.Ny])
                try
                    for slice = 1:nr_slices
                        arr_out(:,:,slice) = process_handle(arr_in(:,:,slice));
                    end
                catch err
                    error(err);
                end
                warning('Array size does not match detector size and was automatically resized.');
            else
                arr_out = ref;
            end
        end
        
        
        function roi = around_y(obj,pixel)
            % AROUND_Y  Description.
            %
            %   roi = around_y(pixel)
            %
            % The following arguments are supported:
            %
            %       pixel:: []
            %           Description.
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   roi:: Description.
            %
            
            roi = round([-pixel pixel] + obj.pby);
        end
        
        function roi = around_z(obj,pixel)
            % AROUND_Z  Description.
            %
            %   roi = around_z(pixel)
            %
            % The following arguments are supported:
            %
            %       pixel:: []
            %           Description.
            %
            % Example:
            %   Example missing.
            %
            % Output arguments:
            %   roi:: Description.
            %
            
            roi = round([-pixel pixel] + obj.pbz);
        end
    end
end
