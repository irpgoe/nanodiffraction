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
classdef spec<handle
    properties
        prepath;
        c;
    end
    methods
        function obj = spec(config)
            obj.c = config;
            obj.prepath = config.path;
        end
        
        function [data,header,scaninfo] = read(obj,specNr,varargin)
            % READ_SPECFILE  Is a handwritten implementation to read 
            % specfiles that are generated at id13. Do not trust that this 
            % function will always work, it has only been used at one 
            % specific beamtime.
            %   
            %   [DATA,HEADER] = READ_SPECFILE(SPECNR, VARARGIN) 
            %
            %   The following options are supported:
            %
            %     specNr:: []
            %       Number of scan in specfile (denoted by #S).
            %
            %     specfile:: [] (not required)
            %       Spec filename.

            switch obj.beamline
                case 'id13'
                    [data,header,scaninfo] = obj.read_id13(specNr,varargin{:});
                case 'p10'
                    [data,header,scaninfo] = obj.read_p10(specNr,varargin{:});
            end
        end
        
        
        function [data,header,scaninfo] = read_id13(obj,specNr,varargin)
            % READ_SPECFILE  Is a handwritten implementation to read 
            % specfiles that are generated at id13. Do not trust that this 
            % function will always work, it has only been used at one 
            % specific beamtime.
            %   
            %   [DATA,HEADER] = READ_SPECFILE(SPECNR, VARARGIN) 
            %
            %   The following options are supported:
            %
            %     specNr:: []
            %       Number of scan in specfile (denoted by #S).
            %
            %     specfile:: [] (not required)
            %       Spec filename.
            
            if nargin > 2
                specfile = varargin{1};
            else 
                specfile = obj.prepath;
            end
            
            % specfile
            fid = fopen(specfile);
            if fid < 0
                fprintf(1,'Debug information \n');
                fprintf(1,'Specfile: %s\n',specfile);
                error('Please choose a valid filename');
            end
            
            % Read until #S is encountered
            tline = fgetl(fid);
            cond = 1;
            while cond < 100
                % split up
                tmp = strsplit(tline,' ');
                if strcmp(tmp(1),'#S')
                    % Compare scan number
                    if strcmp(tmp(2),num2str(specNr))
                        scaninfo = tmp;
                        cond = 100;
                    end
                end
                tline = fgetl(fid);
                cond = cond+1;
            end
            
            % If scan number correct, find #L line
            cond = 1;
            while cond < 100
                % split up
                tmp = strsplit(tline,' ');
                if strcmp(tmp(1),'#L')
                    cond = 100;
                    header = tmp;
                    % make some slight adjustments
                    header(1) = [];header(6) = [];
                end
                tline = fgetl(fid);
                cond = cond+1;
            end
            
            while strcmp(tline(1:2),'#C')
                % skip these lines
                tline = fgetl(fid);
            end
            
            % read data
            ii = 1;
%             fprintf(1,'Reading line %5.0f', 1);
            while ~strcmp(tline,'') && ~strcmp(tline(1:2),'#C') && ii < 100000
%                 fprintf(1,'\b\b\b\b\b%5.0f', ii);
                data(ii,:) = str2double(strsplit(tline));     
                tline = fgetl(fid);
                ii = ii + 1;
            end
%             fprintf(1,'\n');

            fclose(fid);
            
        end
        
        
        
        
        function [data,header,scaninfo] = read_p10(obj,specNr,varargin)
            % READ_SPECFILE  Is a handwritten implementation to read 
            % specfiles that are generated at id13. Do not trust that this 
            % function will always work, it has only been used at one 
            % specific beamtime.
            %   
            %   [DATA,HEADER] = READ_SPECFILE(SPECNR, VARARGIN) 
            %
            %   The following options are supported:
            %
            %     specNr:: []
            %       Number of scan in specfile (denoted by #S).
            %
            %     specfile:: [] (not required)
            %       Spec filename.
            
            if nargin > 2
                specfile = varargin{1};
            else 
                specfile = obj.prepath;
            end

            % specfile
            fid = fopen(specfile);
            if fid < 0
                fprintf(1,'Debug information \n');
                fprintf(1,'Specfile: %s\n',specfile);
                error('Please choose a valid filename');
            end
            
            % Read until #S is encountered
            tline = fgetl(fid);
            cond = 1;
            while cond < 10000
                % split up
                tmp = strsplit(tline,' ');
                if strcmp(tmp(1),'#S')
                    % Compare scan number
                    if strcmp(tmp(2),num2str(specNr))
                        cond = 10000;
                        scaninfo = tmp;
                    end
                end
                tline = fgetl(fid);
                cond = cond+1;
            end
            
            % If scan number correct, find #L line
            cond = 1;
            while cond < 100
                % split up
                tmp = strsplit(tline,' ');
                if strcmp(tmp(1),'#L')
                    cond = 100;
                    header = tmp;
                    % make some slight adjustments
                    header(1) = [];
                end
                tline = fgetl(fid);
                cond = cond+1;
            end
                        
            % read data
            ii = 1;
            fprintf(1,'Reading line %5.0f', 1);
            while ~strcmp(tline,'') && ii < 100000 
                fprintf(1,'\b\b\b\b\b%5.0f', ii);
                data(ii,:) = str2double(strsplit(tline));     
                tline = fgetl(fid);
                if tline < 0
                    ii = 100000;
                end
                
                if size(tline,1) > 1
                    if strcmp(tline(1:2),'#C')
                        ii = 100000;
                    end
                end
                ii = ii + 1;
            end
            fprintf(1,'\n');

            fclose(fid);
            
        end
    end
end