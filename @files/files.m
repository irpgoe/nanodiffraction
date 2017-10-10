classdef files<handle
    % FILES  Class to process file requests for nanodiffraction data.
    % files accepts name value pairs and requires in general the path to
    % certain scan parameters such as the directory of data storage, file
    % naming conventions (prefixes and suffixes) and the type of detector
    % used. Once initialized properly, it will be attached to a
    % nanodiffraction class so that file reading will be done
    % automatically based on the scan frame number.
    %   
    %   DATA = FILES(VARARGIN) creates an files object, with a list of 
    %   arguments provided in VARARGIN. VARARGIN consists of NAME - VALUE 
    %   pairs. NAME is a string describing the option and VALUE can be any 
    %   input. Typically, a specific input is requested by default and an 
    %   error is thrown if the types do not agree.
    %
    %   The following options are supported:
    %
    %   prepath:: [homegroups/Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1] 
    %     path to folder where the data is stored. The prepath depends 
    %     wether you choose p10 or id13 as beamline. 
    %       id13: it is the path right before the data directory, hence
    %       prepath/data_???_master.h5, usually ending with
    %       .../DATA/AUTO-TRANSFER/eiger1
    %       p10: it is the path before the /detectors/ directory, hence
    %       prepath/detectors/<detector>/<newfile>/data_???_master.h5
    %
    %   beamline:: [id13]
    %     Beamline, where data was recorded. Currently, only 'id13' and
    %     'p10' are supported.
    %
    %   newfile:: [detdistcalib] 
    %     name (prefix) of scan(s)
    %
    %   detector:: [eiger] 
    %     Type of detector used in the experiment.
    %     currently, only 'eiger' and 'pilatus' are supported
    %
    %   specNr:: [1]
    %      Number of the scan as marked in the spec file.
    %
    %   fluoNr:: [1]
    %     suffix of fluorescence data usually stored in .dat files
    %
    %   fn:: [1]
    %     file (frame) number
    %
    %   fluoNewfile:: []
    %     the fluorescence detector might have its own file naming 
    %     convention
    %
    %   fluoDetector:: [xia]
    %     currently only 'xia' (ID13) is supported
    %
    %   fluoFramesPerFile:: [101]
    %     usually the number of scan points along the fast axis are given
    %
    %   eigerNr:: [210]
    %     Number of the master file (usually ending with _master.h5)
    %
    %   Example::
    %     fp = files( 'beamline','id13',...
    %       'prepath',[prepath 'Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1'],...
    %       'newfile','herz2_roi2',...
    %       'detector','eiger',...
    %       'eigerNr',201);
    %
    %   Properties::
    %      <data> - p (holds all the necessary parameters)
    %             - specpath (path to spec files)
    %             - fluopath (path to fluorescence data)
    %             - debug (debug level) [0]
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
        p = struct('fnConvention',5,'fileFormat','cbf'); % parameter definitions
        specpath;
        fluopath;
        detectors = struct();
        debug = 0;
        
%         % eigerspecific file reading, 
%         % expects str1 : _master.h5 file
%         %         str2 : /entry/data/data_%f6.0 file number
%         eiger_read_total = @(str1,str2) double(h5read( str1, str2 ));
%         eiger_read_single = @(str1,str2,vec1,vec2) double(h5read( str1, str2, vec1, vec2 ));
            
    end
    methods

        function obj = files(varargin)
            
        	p = inputParser;
            % set defaults and add optional arguments
            addOptional(p,'prepath','homegroups/Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1');
            addOptional(p,'newfile','detdistcalib');
            addOptional(p,'detector','eiger');
            addOptional(p,'beamline','id13');
            addOptional(p,'eigerNr',210);
            addOptional(p,'specNr',1);
            addOptional(p,'fluoNr',1);
            addOptional(p,'fn',1);
            addOptional(p,'fluoNewfile','');
            addOptional(p,'fluoDetector','xia');
            addOptional(p,'fluoFramesPerFile',101);
            addOptional(p,'framesPerFile',1);
            addOptional(p,'fluopath','');
            addOptional(p,'specpath','');
            
            % parse arguments
            parse(p,varargin{:});
            % copy p.Results to obj.p
            names = fieldnames(p.Results);
            for ii = 1:numel(names)
                obj.p.(names{ii}) = p.Results.(names{ii});
            end
            
            % check path
            [obj.p.prepath,obj.p.slash] = obj.which_slash(obj.p.prepath);
            obj.p.prepath = obj.remove_leading_trailing_slashes(obj.p.prepath);
            
            
            
            if contains(p.Results.detector,'eiger')
                obj.detectors.eiger = eiger(obj.p.beamline,obj.p.prepath,obj.p.newfile,obj.p.eigerNr,obj.p.framesPerFile,obj.p.slash);
%                 obj.set_eiger_scan(obj.p.newfile,obj.p.eigerNr);
                obj.repair_frames_per_file();
            end
            
            if contains(p.Results.detector,'rayonix')
                obj.detectors.rayonix = rayonix(obj.p.beamline,obj.p.prepath,obj.p.newfile,obj.p.slash);            
            end
            
            if contains(p.Results.detector,'pilatus')
                obj.detectors.pilatus = pilatus(obj.p.beamline,obj.p.prepath,obj.p.newfile,obj.p.slash);            
            end
            
            if contains(p.Results.detector,'xia')
                obj.detectors.xia = xia(obj.p.beamline,obj.p.fluopath,obj.p.fluoNewfile,obj.p.fluoNr,obj.p.fluoFramesPerFile,obj.p.slash);
            end            
            
            if contains(p.Results.detector,'spec')
                obj.detectors.spec = spec(obj.p.beamline,obj.p.prepath,obj.p.newfile,obj.p.slash);
            end
        end
        
        
        function set_frames_per_file(obj,fpf)
            % SET_FRAMES_PER_FILE  sets the number of frames per file if
            % data is stored in a container format such as hdf5.
            %   
            %   SET_FRAMES_PER_FILE(FPF) 
            %
            %   The following options are supported:
            %
            %     fpf:: []
            %       Frames per file.
            obj.detectors.eiger.framesPerFile = fpf;
        end
        
        
        function set_eiger_scan(obj,newfile,eigerNr)
            % SET_EIGER_SCAN  updates the eiger module,
            %   and automatically determines the frames per file.
            %   
            %   SET_EIGER_SCAN(newfile, eigerNr) 
            %
            %   The following options are supported:
            %
            %     newfile:: []
            %       Scan name.
            %
            %     eigerNr:: []
            %       Scan number.
            obj.p.newfile = newfile;
            obj.p.eigerNr = eigerNr;
            obj.repair_frames_per_file();
            obj.detectors.eiger = eiger(obj.p.beamline,obj.p.prepath,newfile,eigerNr,obj.p.framesPerFile,obj.p.slash);
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
            obj.p.newfile = newfile;
            obj.detectors.pilatus = pilatus(obj.p.beamline,obj.p.prepath,newfile,obj.p.slash);
        end
        
        
        function set_xia_scan(obj,fluoNewfile,fluoNr,fluoFramesPerFile)
            % SET_XIA_SCAN  updates the xia module.
            %   
            %   SET_XIA_SCAN(newfile, fluoNr, fpf) 
            %
            %   The following options are supported:
            %
            %     newfile:: []
            %       Scan name.
            %
            %     fluoNr:: []
            %       Scan number.
            %
            %     fpf:: []
            %       Scan number.
            obj.p.fluoNewfile = fluoNewfile;
            obj.p.fluoNr = fluoNr;
            obj.p.fluoFramesPerFile = fluoFramesPerFile;
            obj.detectors.eiger = eiger(obj.p.beamline,obj.p.fluopath,fluoNewfile,fluoNr,fluoFramesPerFile,obj.p.slash);
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
            obj.p.newfile = newfile;
            obj.detectors.spec = spec(obj.p.beamline,obj.p.prepath,newfile,obj.p.slash);
        end
        

        function [out_string] = remove_leading_trailing_slashes(obj,string)
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
                
        
        function repair_frames_per_file(obj)
            % REPAIR_FRAMES_PER_FILE  calculates frames per file when using 
            % the h5 container format e.g. with the eiger detector
            %
            %   REPAIR_FRAMES_PER_FILE()
            
            % starting point
            n = 1;
            m = 100;
            fpf = 1;
            cond = 1;
            while cond
                
                % try to read
                try double(h5read(obj.compile_path{1},...
                        ['/entry/data/data_' sprintf('%06i',1)],...
                        [1 1 n],...
                        [2070 2167 1]));
                catch
                    % we have gone out of bounds
                    n = n - m; % go back to value before
                    m = m / 10; % reduce sampling rate
                    if m < 1 
                        cond = 0;
                    end
                end
                % last working condition
                fpf = n;
                % we are still within bounds, hence step forward
                n = n + m;
                
            end
            obj.p.framesPerFile = fpf;
            obj.set_frames_per_file(fpf);
            
            % debugging
            % fprintf(1,'Frames per file and file number: %d, %d\n',obj.p.framesPerFile,obj.p.fn);
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
            %       Frame number. If frame number is not given, obj.p.fn
            %       will be used
            
            
            % if filenumber is provided, use it
            if nargin >= 2
                obj.p.fn = varargin{1};
            end
                
            % remove leading or trailing slashes
            if ~isempty(obj.p.prepath) && ~isempty(obj.p.newfile) && ~isempty(obj.p.detector) && obj.p.fn >= 0
                if strcmp(obj.p.prepath(1),'/') || strcmp(obj.p.prepath(1),'\')
                    obj.p.prepath(1) = []; 
                elseif strcmp(obj.p.prepath(end),'/') || strcmp(obj.p.prepath(end),'\')
                    obj.p.prepath(end) = []; 
                end
            else
                % return error
                disp('Prepath, Newfile or detector might not be set. Please set these values accordingly before continuing!');
                return
            end
            
            % slash
            [obj.p.prepath, sl] = obj.which_slash(obj.p.prepath);            
            
            % beamline-specific naming conventions
            switch obj.p.beamline
                case 'p10'
                    switch obj.p.detector
                        case 'eiger'
                            compiledPath{1} = [sl obj.p.prepath sl 'detectors' sl 'eiger' sl obj.p.newfile sl obj.p.newfile '_' sprintf('%05i',obj.p.eigerNr) '_' 'master.h5'];
                        case 'pilatus'
                            compiledPath{1} = [sl obj.p.prepath sl 'detectors' sl 'p300' sl obj.p.newfile sl obj.p.newfile '_' sprintf('%05i',obj.p.fn) '.cbf'];
                    end
                case 'id13'
                    switch obj.p.detector
                        case 'eiger'
                            compiledPath{1} = [sl obj.p.prepath sl obj.p.newfile '_' num2str(obj.p.eigerNr) '_' 'master.h5'];
                    end
                case 'id02'
                    compiledPath{1} = [sl obj.p.prepath sl obj.p.newfile '_' sprintf('%05i',obj.p.fn) '_raw.edf'];
                case 'none'
                    switch obj.p.detector
                        case 'eiger'
%                             compiledPath{1} = [sl obj.p.prepath sl obj.p.newfile '_' num2str(obj.p.eigerNr) '_' 'master.h5'];
                            compiledPath{1} = [sl obj.p.prepath sl obj.p.newfile '_' sprintf('%05i',obj.p.eigerNr) '_' 'master.h5'];
                    end
            end
            
            % detector-specific
            if strcmp(obj.p.detector,'eiger')
                index = floor(double(obj.p.fn-1)/double(obj.p.framesPerFile))+1;
                remainder = mod(double(obj.p.fn-1),double(obj.p.framesPerFile))+1;
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
                data = obj.detectors.(obj.p.detector).read(varargin{:});
            else
                % otherwise, use last frame number that was in use.
                data = obj.detectors.(obj.p.detector).read(obj.p.fn);
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
            if ~isempty(obj.p.prepath) && ~isempty(obj.p.newfile) && ~isempty(obj.p.detector) && obj.p.fn >= 0
                if strcmp(obj.p.prepath(1),'/') || strcmp(obj.p.prepath(1),'\')
                    obj.p.prepath(1) = []; 
                elseif strcmp(obj.p.prepath(end),'/') || strcmp(obj.p.prepath(end),'\')
                    obj.p.prepath(end) = []; 
                end
                
            else
                % return error
                disp('Prepath, Newfile or detector might not be set. Please set these values accordingly before continuing!');
                return
            end
            
            % slash
            [obj.p.prepath, sl] = obj.which_slash(obj.p.prepath);            
            
            % slice variables
            fpf = obj.p.framesPerFile;
            fileFormat = obj.p.fileFormat;
            
            % then compile function
            switch obj.p.beamline
                case 'p10'
                    switch obj.p.detector
                        case 'eiger'
                            staticPath = [sl obj.p.prepath sl 'detectors' sl 'eiger' sl obj.p.newfile sl obj.p.newfile '_' sprintf('%05i',obj.p.eigerNr) '_master.h5'];
                            
                            f = @(fn) flipud(rot90(double(h5read(staticPath,...
                                ['/entry/data/data_' sprintf('%06i',floor(double(fn-1)/double(fpf))+1)],...
                                [starty startx double(mod(double(fn-1),double(fpf))+1)],...
                                [lengthy lengthx 1]))));
                        case 'pilatus'
                            staticPath = [sl obj.p.prepath sl 'detectors' sl 'p300' sl obj.p.newfile sl obj.p.newfile '_'];
                            f = @(fn) flipud(imagereader_P10([staticPath sprintf('%05i',fn) '.cbf' ],'p300'));
                    end
                case 'id13'
                    switch obj.p.detector
                        case 'eiger'
                            staticPath = [sl obj.p.prepath sl obj.p.newfile '_' num2str(obj.p.eigerNr) '_master.h5'];

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
        
        
%         function result = read_from_dada(~, varargin)
%         % READ_FROM_DATA  reads (preprocessed) data from dada. 
%         % 
%         %   DATA = READ_FROM_DADA(VARARGIN)
%         %
%         %   The following arguments are accepted:
%         %
%         %     INSTRUMENT:: [ESRF]
%         %       As in dada.
%         %
%         %     EXPERIMENT:: [ls2522]
%         %       As in dada.
%         %
%         %     DETECTOR:: [eiger1]
%         %       As in dada.
%         %
%         %     SAMPLE:: [n05_ms532_uu09_roi6new]
%         %       As in dada.
%         %
%         %     FILENUMBER:: [11]
%         %       As in dada.
%         %
%         %     FRAMENUMBER:: [1]
%         %       As in dada.
%         %
%         % In addition you could specify a ROI or BINNING
%         %     ROI:: []
%         %       e.g. '1090,1111,200,200'
%         % 	  BINNING:: ['1,1']
%         %       e.g. '1,1'
%         
%             
%             % parse input
%             pin = inputParser;
%             addOptional(pin,'experiment','ls2522');
%             addOptional(pin,'detector','eiger1');
%             addOptional(pin,'instrument','ESRF');
%             addOptional(pin,'sample','n05_ms532_uu09_roi6new');
%             addOptional(pin,'filenumber',11);
%             addOptional(pin,'framenumber',1); 
%             addOptional(pin,'roi',''); 
%             addOptional(pin,'binning','1,1');
%             parse(pin,varargin{:});
%             
%             if (nargin == 1)
%                disp('Showing you an example scan');
%                disp(pin.Results)
%             end
%             
%             % rename
%             experiment = pin.Results.experiment;
%             detector = pin.Results.detector;
%             instrument = pin.Results.instrument;
%             sample = pin.Results.sample;
%             filenumber = pin.Results.filenumber;
%             framenumber = pin.Results.framenumber;
%             
%             % preprocessing
%             roi = pin.Results.roi;
%             binning = pin.Results.binning;
%             
%             % get data from data
%             tic;
%             dada = 'http://dada-devel.roentgen.physik.uni-goettingen.de';
%             url  =  [dada '/show/' instrument '/' experiment '/' detector '/' sample '/' num2str(filenumber) '/' num2str(framenumber) ['/?format=double&roi=' roi '&binning=' binning]];
%             [data,extra] = urlread2(url);
%             dims = str2num(extra.allHeaders.X_dimensions{1});
%             data = typecast(cast(data, 'uint8'), 'double');
%             data = reshape(data, dims);
%             toc;
% 
%             d = zeros(dims)-1;
%             d(data>0) = (data(data>0));
% %             d(data>0) = log10(data(data>0));
% %             imagesc(flipud(rot90(d))); caxis([-1 3]); colormap gray; colorbar;
% %             fprintf('compare to http://dada-devel/show/GINIX/run55/eiger/simone/17/1?scale=log&palette=bw&binning=4,4&cmin=0&cmax=1000\n');
% %             toc;
% 
%             % final signal
%             result = flipud(rot90(d));
% 
%         end
        
                
        function firstFile = set_first_file(obj,varargin)
            % SET_FIRST_FILE  sets the scan information for all subsequent 
            % analyses. For p10, one could use spec_get_scan_info. This
            % function is necessarily to be called when e.g. images were
            % acquired in eigerslow mode, because then single images are
            % placed into a single big h5 file.
            %
            %    FIRSTFILE = SET_FIRST_FILE(VARARGIN)
            %
            %    The following arguments are accepted:
            %
            %      beamline:: [p10]
            %        The beamline in use.
            %
            %      scanNo:: [1]
            %        The number of the scan.
            %
            
            if nargin == 2
                obj.p.firstFile = varargin{1};
            else
            
                % parse input
                pin = inputParser;
                addOptional(pin,'beamline','p10');
                addOptional(pin,'scanNo',1);
                parse(pin,varargin{:});

                switch pin.Results.beamline
                    case 'p10'
                        % scan details
                        scan_details = spec_get_scan_info(pin.Results.scanNo, obj.p.newfile,obj.specpath);

                        % file numbers for each scanpoint
                        firstFile  =  scan_details.first_file;
                        obj.p.firstFile = firstFile;

                    case 'id13'
                        % file numbers for each scanpoint
                        firstFile = 1;

                        % save to object
                        obj.p.firstFile = firstFile;                       

                end
            end
        end
    end
end