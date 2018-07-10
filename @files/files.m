classdef files<handle
    % FILES  Class to process file requests for nanodiffraction data.
    % FILES accepts name value pairs and requires in general the path to
    % certain scan parameters such as the directory of data storage, file
    % naming conventions (prefixes and suffixes) and the type of detector
    % used. Once initialized properly, it should be linked to a
    % nanodiffraction class so that file reading will be done
    % automatically based on the scan frame number. Use link() to connect a
    % FILES class to a NANODIFFRACTION class.
    %   
    %   DATA = FILES(VARARGIN) creates a files object, with a list of 
    %   arguments provided in VARARGIN. VARARGIN consists of NAME - VALUE 
    %   pairs. NAME is a string describing the option and VALUE can be any 
    %   input. Typically, a specific input is requested by default and an 
    %   error is thrown if the types do not agree.
    %
    %   The following options are supported:
    %
    %   prepath:: [homegroups/Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1] 
    %     Path to folder where the data is stored. The prepath depends 
    %     on the beamline you choose. Check help path_by_beamline for more
    %     information on the default settings.
    %
    %   beamline:: [id13]
    %     The beamline argument is necessary, because it is referred to
    %     when a certain standard file naming convention should be used.
    %     See the help on path_by_beamline for more information.
    %
    %       Options: 'id13'|'p10'|'id02'|'bessy'|'none' 
    %
    %   newfile:: [detdistcalib] 
    %     File name prefix.
    %
    %   detector:: [eiger] 
    %     Type of detector used in the experiment.
    %     Each detector is currently contained in a separate module. 
    %
    %       Options: 'eiger'|'pilatus'|'rayonix'|'xia'
    %
    %   fn:: [1]
    %     File number. Note, that some detectors start the numbering scheme
    %     with 1, while others start with 0. 
    %
    %   fpf:: [101]
    %     Frames per file. If a container format is used that collects a
    %     given number (the fpf) of scattering patterns, this number of
    %     frames per container is termed "frames per file". E.g. the Eiger
    %     detector uses the hdf5 container format. fpf can be an arbitrary
    %     number, but usually it is the number of scan points along the 
    %     fast axis.
    %
    %   scan:: [210]
    %     Number of the scan. Be aware, that SPEC or other logging software
    %     might use a different numbering scheme than the detector. The
    %     scan number always refers to the number that is appended to the
    %     filename.
    %
    %   Example::
    %     fp = files( 'beamline','id13',...
    %       'prepath',[prepath 'Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1'],...
    %       'newfile','herz2_roi2',...
    %       'detector','eiger',...
    %       'scan',201);
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
        
        prepath;    % see help files
        newfile;    % see help files
        beamline;   % see help files
        detector;   % see help files
        scan;       % see help files
        fn;         % see help files
        fpf;        % see help files
        firstFile;  % file number of the first frame in a scan
        slash;      % either '/' or '\' based on the system in use
        
        % internal use only
        detectors = struct();
        debug = 0;  
            
    end
    methods

        function obj = files(varargin)
            
        	p = inputParser;
            % set defaults and add optional arguments
            addOptional(p,'prepath','homegroups/Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1');
            addOptional(p,'newfile','detdistcalib');
            addOptional(p,'detector','eiger');
            addOptional(p,'beamline','id13');
            addOptional(p,'scan',210);
            addOptional(p,'fn',1);
            addOptional(p,'fpf',101);
            
            % parse arguments
            parse(p,varargin{:});
            % copy p.Results to obj.p
            names = fieldnames(p.Results);
            for ii = 1:numel(names)
                obj.(names{ii}) = p.Results.(names{ii});
            end
            
            
            % check path
            [obj.prepath,obj.slash] = obj.which_slash(obj.prepath);
            obj.prepath = obj.remove_leading_trailing_slashes(obj.prepath);

            
            % get file information by beamline and detector
            config = obj.path_by_beamline(obj.beamline,p.Results.detector,obj.prepath,obj.newfile,obj.scan);

            switch p.Results.detector
                case 'eiger'
                    obj.detectors.eiger = eiger(config,obj.fpf);
                case 'rayonix'
                    obj.detectors.rayonix = rayonix(config);            
                case 'pilatus'
                    obj.detectors.pilatus = pilatus(config);            
                case 'pirra'
                    obj.detectors.pirra = pirra(config);         
                case 'talos'
                    obj.detectors.talos = talos(config);                       
                case 'xia'
                    obj.detectors.xia = xia(config,obj.fpf);
                case 'spec'
                    obj.detectors.spec = spec(config);
            end
        end
        
        
        function [config] = path_by_beamline(obj,beamline,detector,prepath,newfile,scannr)
            % PATH_BY_BEAMLINE  contains the appropriate folder structure
            % for each beamline. This file can be adapted to the specific
            % needs. 
            %   
            %   CONFIGURATION = PATH_BY_BEAMLINE(BEAMLINE,DETECTOR,PREPATH,NEWFILE,SCANNR) 
            %
            %   The following options are supported:
            %
            %     beamline:: [no default] 
            %       Beamline short name. As described in help files.
            %
            %       Options: 'p10'|'id13'|'id02'|'bessy'
            %
            %     detector: [no default] 
            %       Detector. As described in help files.
            %
            %       Options:
            %       'pilatus'|'eiger'|'pirra'|'rayonix'|'xia'|'spec'
            %
            %     prepath: [no default]
            %       Path to data directory. As described in help files.
            %
            %     newfile: [no default]
            %       Name of scan. As described in help files.
            %
            %     scannr: [no default]
            %       Scan number. As described in help files.
            %
    
            % shorthand
            sl = obj.slash;

            setting = [detector '_' beamline];
            config = struct();
            switch setting
                case 'eiger_p10'
                    config.path = [sl prepath...
                        sl 'detectors'...
                        sl 'eiger'...
                        sl newfile...
                        sl newfile '_'...
                        sprintf('%05i',scannr) '_master.h5'];
                case 'eiger_id13'
                    config.path = [sl prepath...
                        sl newfile '_'...
                        num2str(scannr) '_master.h5'];
                case 'eiger_none'
                    config.path = [sl prepath...
                        sl newfile '_'...
                        sprintf('%i',scannr) '_master.h5'];
                case 'pilatus_p10'
                    if ~ispc
                        prepath = [sl prepath];
                    end
                    config.path = [prepath...
                        sl 'detectors'...
                        sl 'p300'...
                        sl newfile...
                        sl newfile '_'];
                    config.filename = @(pre,n) [pre sprintf('%05i',n) '.cbf'];
                case 'pirra_p10'
                    if ~ispc
                        prepath = [sl prepath];
                    end
                    config.path = [prepath...
                        sl 'detectors'...
                        sl 'pirra'...
                        sl newfile...
                        sl newfile '_'];
                    config.filename = @(pre,n) [pre sprintf('%i',n) '.tif'];    
                case 'talos_none'
                    if ~ispc
                        prepath = [sl prepath];
                    end                    
                    config.path = [prepath...
                        sl 'detectors'...
                        sl 'talos'...
                        sl newfile...
                        sl newfile '_'];
                    config.filename = @(pre,n) [pre sprintf('%04i',n) '.tif']; 
                case 'rayonix_id02'
                    config.path = [sl prepath...
                        sl newfile '_'];        
                    config.filename = @(pre,n) [pre sprintf('%05i',n) '_raw.edf'];
                    config.N = 960;
                case 'rayonix_bessy'
                    config.path = [sl prepath...
                        sl newfile '_' sprintf('%05i',scannr) '_'];        
                    config.filename = @(pre,n) [pre sprintf('%05i',n) '.mccd'];
                    config.N = 3072;
                case 'spec_id13'        
                    config.path = [sl prepath...
                        sl newfile...
                        sl newfile '.dat'];     
                    config.beamline = 'id13';
                case 'spec_p10'
                    config.path = [sl prepath...
                        sl newfile];   
                    config.beamline = 'p10';
                case 'xia_id13'        
                    config.path = [sl prepath...
                        sl newfile...
                        sl 'xia'...
                        sl 'xia'...
                        '_xia02_' sprintf('%04i',scannr) '_0000_'];     
                    config.filename = @(pre,n) [pre sprintf('%04i',n) '.edf'];
                case 'default'
                    config.path = [sl prepath...
                        sl newfile '_'...
                        sprintf('%05i',scannr) '_master.h5'];                    
            end            
        end          
        
        
        function set_eiger_scan(obj,newfile,scan)
            % SET_EIGER_SCAN  updates the eiger module.
            %   
            %   SET_EIGER_SCAN(newfile, scan) 
            %
            %   The following options are supported:
            %
            %     newfile:: []
            %       Scan name.
            %
            %     scan:: []
            %       Scan number.
            
            obj.newfile = newfile;
            obj.scan = scan;
            
            config = obj.path_by_beamline(obj.beamline,'eiger',obj.prepath,obj.newfile,obj.scan);
            obj.detectors.eiger = eiger(config,obj.fpf);
        end
        
        
        function set_pilatus_scan(obj,newfile)
            % SET_PILATUS_SCAN  updates the pilatus module.
            %   
            %   SET_PILATUS_SCAN(newfile) 
            %
            %   The following options are supported:
            %
            %     newfile:: []
            %       Scan name.
            
            obj.newfile = newfile;
            
            config = obj.path_by_beamline(obj.beamline,'pilatus',obj.prepath,obj.newfile,obj.scan);
            obj.detectors.pilatus = pilatus(config);
        end
        
        
        function set_xia_scan(obj,newfile,scan,fpf)
            % SET_XIA_SCAN  updates the xia module.
            %   
            %   SET_XIA_SCAN(newfile, scan, fpf) 
            %
            %   The following options are supported:
            %
            %     newfile:: []
            %       Scan name.
            %
            %     scan:: []
            %       Scan number.
            %
            %     fpf:: []
            %       Scan number.
            
            obj.newfile = newfile;
            obj.scan = scan;
            obj.fpf = fpf;
            
            config = obj.path_by_beamline(obj.beamline,'xia',obj.prepath,obj.newfile,obj.scan);
            obj.detectors.xia = xia(config,fpf);
        end
        
        
        function set_spec_scan(obj,newfile)
            % SET_SPEC_SCAN  updates the spec module.
            %   
            %   SET_SPEC_SCAN(newfile) 
            %
            %   The following options are supported:
            %
            %     newfile:: []
            %       Scan name.
            
            obj.newfile = newfile;
            
            config = obj.path_by_beamline(obj.beamline,'spec',obj.prepath,obj.newfile,obj.scan);
            obj.detectors.spec = spec(config);
        end
        

        function [out_string] = remove_leading_trailing_slashes(~,string)
            % REMOVE_LEADING_TRAILING_SLASHES  removes any leading or
            % trailing slashes from a given string.
            %   
            %   REMOVE_LEADING_TRAILING_SLASHES(STRING) 
            %
            %   The following options are supported:
            %
            %     string:: []
            %       Input string.
            
            if isempty(string)
                out_string = string;
                return;
            end
            
            if strcmp(string(1),'/') || strcmp(string(1),'\')
                string(1) = []; 
            end
            if strcmp(string(end),'/') || strcmp(string(end),'\')
                string(end) = []; 
            else
                % nothing to do
            end
            out_string = string;
        end
        
        
        function [compiledPath] = compile_path(obj,varargin)
            % COMPILE_PATH  generates a complete path or list of paths to a
            % file. This function is called whenever a file (frame) should 
            % be read. It can provide the option to generate the path to a
            % single file or multiple files based on the filenumber
            % provided. 
            %   
            %   [COMPILEDPATH] = COMPILE_PATH(VARARGIN) 
            %
            %   The following options are supported:
            %
            %     Argument 1::
            %       Frame number. If frame number is not given, obj.fn
            %       will be used
            
            
            % if filenumber is provided, use it
            if nargin >= 2
                obj.fn = varargin{1};
            end
                
            % remove leading or trailing slashes
            if ~isempty(obj.prepath) && ~isempty(obj.newfile) && ~isempty(obj.detector) && obj.fn >= 0
                if strcmp(obj.prepath(1),'/') || strcmp(obj.prepath(1),'\')
                    obj.prepath(1) = []; 
                elseif strcmp(obj.prepath(end),'/') || strcmp(obj.prepath(end),'\')
                    obj.prepath(end) = []; 
                end
            else
                % return error
                disp('Prepath, Newfile or detector might not be set. Please set these values accordingly before continuing!');
                return
            end
            
            % slash
            [obj.prepath, sl] = obj.which_slash(obj.prepath);            
            
            % beamline-specific naming conventions
            switch obj.beamline
                case 'p10'
                    switch obj.detector
                        case 'eiger'
                            compiledPath{1} = [sl obj.prepath sl 'detectors' sl 'eiger' sl obj.newfile sl obj.newfile '_' sprintf('%05i',obj.scan) '_' 'master.h5'];
                        case 'pilatus'
                            compiledPath{1} = [sl obj.prepath sl 'detectors' sl 'p300' sl obj.newfile sl obj.newfile '_' sprintf('%05i',obj.fn) '.cbf'];
                    end
                case 'id13'
                    switch obj.detector
                        case 'eiger'
                            compiledPath{1} = [sl obj.prepath sl obj.newfile '_' num2str(obj.scan) '_' 'master.h5'];
                    end
                case 'id02'
                    compiledPath{1} = [sl obj.prepath sl obj.newfile '_' sprintf('%05i',obj.fn) '_raw.edf'];
                case 'none'
                    switch obj.detector
                        case 'eiger'
%                             compiledPath{1} = [sl obj.prepath sl obj.newfile '_' num2str(obj.scan) '_' 'master.h5'];
                            compiledPath{1} = [sl obj.prepath sl obj.newfile '_' sprintf('%05i',obj.scan) '_' 'master.h5'];
                    end
            end
            
            % detector-specific
            if strcmp(obj.detector,'eiger')
                index = floor(double(obj.fn-1)/double(obj.fpf))+1;
                remainder = mod(double(obj.fn-1),double(obj.fpf))+1;
                compiledPath{2} = ['/entry/data/data_' sprintf('%06i',index)];
                compiledPath{3} = double(remainder);
            end
                
        end

        function data = read(obj,varargin)
            % READ  reads data frame from the requested detector module.
            %   
            %   READ(fn) 
            %
            %   The following options are supported:
            %
            %     fn:: []
            %       Frame number.
            
            % if frame number is given
            if nargin > 1
                data = obj.detectors.(obj.detector).read(varargin{:});
            else
                % otherwise, use last frame number that was in use.
                data = obj.detectors.(obj.detector).read(obj.fn);
            end
        end
            
            
        function f = compile_static_function(obj)
            % COMPILE_STATIC_FUNCTION  compiles a function handle with
            % static parameters that can be used in a parallelization task
            % to distribute a simple function handle to the workers
            %
            %   FUNCTIONHANDLE = COMPILE_STATIC_FUNCTION() creates a
            %   function handle based on the current setup of the files
            %   object.
            %
            
            % by default
            starty = 1;
            lengthy = 2070;
            startx = 1;
            lengthx = 2167;
                
            % remove leading or trailing slashes
            if ~isempty(obj.prepath) && ~isempty(obj.newfile) && ~isempty(obj.detector) && obj.fn >= 0
                if strcmp(obj.prepath(1),'/') || strcmp(obj.prepath(1),'\')
                    obj.prepath(1) = []; 
                elseif strcmp(obj.prepath(end),'/') || strcmp(obj.prepath(end),'\')
                    obj.prepath(end) = []; 
                end
                
            else
                % return error
                disp('Prepath, Newfile or detector might not be set. Please set these values accordingly before continuing!');
                return
            end
            
            % slash
            [obj.prepath, sl] = obj.which_slash(obj.prepath);            
            
            % slice variables
            fpf = obj.fpf;
            
            % then compile function
            switch obj.beamline
                case 'p10'
                    switch obj.detector
                        case 'eiger'
                            staticPath = [sl obj.prepath sl 'detectors' sl 'eiger' sl obj.newfile sl obj.newfile '_' sprintf('%05i',obj.scan) '_master.h5'];
                            
                            f = @(fn) flipud(rot90(double(h5read(staticPath,...
                                ['/entry/data/data_' sprintf('%06i',floor(double(fn-1)/double(fpf))+1)],...
                                [starty startx double(mod(double(fn-1),double(fpf))+1)],...
                                [lengthy lengthx 1]))));
                        case 'pilatus'
                            staticPath = [sl obj.prepath sl 'detectors' sl 'p300' sl obj.newfile sl obj.newfile '_'];
                            f = @(fn) flipud(imagereader_P10([staticPath sprintf('%05i',fn) '.cbf' ],'p300'));
                    end
                case 'id13'
                    switch obj.detector
                        case 'eiger'
                            staticPath = [sl obj.prepath sl obj.newfile '_' num2str(obj.scan) '_master.h5'];

                            f = @(fn) flipud(rot90(double(h5read(staticPath,...
                                ['/entry/data/data_' sprintf('%06i',floor(double(fn-1)/double(fpf))+1)],...
                                [starty startx double(mod(double(fn-1),double(fpf))+1)],...
                                [lengthy lengthx 1]))));
                    end
            end
            
            disp('The following function has been compiled:')
            f
        end
           
        
        function [string, slash] = which_slash(~, string)
            % WHICH_SLASH  checks for slashes in a string. If the string 
            % contains slashes they have to conform to naming convenctions 
            % based on the system (unix, windows) in use
            %
            %   [STRING,SLASH] = WHICH_SLASH(STRING)
            %
            %   The following arguments are accepted:
            %
            %     STRING::
            %        This can be any string.
            
            
            if ispc
                slash = '\';
                string = strrep(string, '/', '\');
            elseif isunix
                slash = '/';
                string = strrep(string, '\', '/');
            elseif ismac
                slash = '/';
                string = strrep(string, '\', '/');
            else
                slash = '/';
            end
        end
        
        
        function result = read_from_dada(~, varargin)
        % READ_FROM_DATA  reads (preprocessed) data from dada. 
        % 
        %   DATA = READ_FROM_DADA(VARARGIN)
        %
        %   The following arguments are accepted:
        %
        %     INSTRUMENT:: [ESRF]
        %       As in dada.
        %
        %     EXPERIMENT:: [ls2522]
        %       As in dada.
        %
        %     DETECTOR:: [eiger1]
        %       As in dada.
        %
        %     SAMPLE:: [n05_ms532_uu09_roi6new]
        %       As in dada.
        %
        %     FILENUMBER:: [11]
        %       As in dada.
        %
        %     FRAMENUMBER:: [1]
        %       As in dada.
        %
        % In addition you could specify a ROI or BINNING
        %     ROI:: []
        %       e.g. '1090,1111,200,200'
        % 	  BINNING:: ['1,1']
        %       e.g. '1,1'
        
            
            % parse input
            pin = inputParser;
            addOptional(pin,'experiment','ls2522');
            addOptional(pin,'detector','eiger1');
            addOptional(pin,'instrument','ESRF');
            addOptional(pin,'sample','n05_ms532_uu09_roi6new');
            addOptional(pin,'filenumber',11);
            addOptional(pin,'framenumber',1); 
            addOptional(pin,'roi',''); 
            addOptional(pin,'binning','1,1');
            parse(pin,varargin{:});
            
            if (nargin == 1)
               disp('Showing you an example scan');
               disp(pin.Results)
            end
            
            % rename
            experiment = pin.Results.experiment;
            detector = pin.Results.detector;
            instrument = pin.Results.instrument;
            sample = pin.Results.sample;
            filenumber = pin.Results.filenumber;
            framenumber = pin.Results.framenumber;
            
            % preprocessing
            roi = pin.Results.roi;
            binning = pin.Results.binning;
            
            % get data from data
            tic;
            dada = 'http://dada-devel.roentgen.physik.uni-goettingen.de';
            url  =  [dada '/show/' instrument '/' experiment '/' detector '/' sample '/' num2str(filenumber) '/' num2str(framenumber) ['/?roi=' roi '&binning=' binning]];
%             url  =  [dada '/show/' instrument '/' experiment '/' detector '/' sample '/' num2str(filenumber) '/' num2str(framenumber) ['/?format=double&roi=' roi '&binning=' binning]];
            [data,extra] = urlread2(url);
            dims = str2num(extra.allHeaders.X_dimensions{1});
            data = typecast(cast(data, 'uint8'), 'double');
            data = reshape(data, dims);
            toc;

            d = zeros(dims)-1;
            d(data>0) = (data(data>0));
%             d(data>0) = log10(data(data>0));
%             imagesc(flipud(rot90(d))); caxis([-1 3]); colormap gray; colorbar;
%             fprintf('compare to http://dada-devel/show/GINIX/run55/eiger/simone/17/1?scale=log&palette=bw&binning=4,4&cmin=0&cmax=1000\n');
%             toc;

            % final signal
            result = flipud(rot90(d));

        end
    end
end