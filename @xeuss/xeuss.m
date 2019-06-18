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
classdef xeuss<handle
    properties
        path;
        size;
        Ny;
        Nz;
        filename;
    end
    methods
        function obj = xeuss(config)
            
            defaults = struct('path','',...
                'size',[981 1043],...
                'filename',@(pre,n) [pre sprintf('%05i',n) '.edf']);
            config = update_defaults(defaults,config);

            % get configuration
            [obj.path, obj.size, obj.filename] = split_struct(config,{'path','size','filename'});
            obj.Ny = obj.size(1);
            obj.Nz = obj.size(2);
        end
        function data = read(obj,fn)
            % compile filepath
            filepath = obj.filename(obj.path,fn);
            
            % open file
            fid = fopen(filepath);
            if fid == -1
                error(['File not found: ' filepath]);
            end
            fseek(fid, 3*512,'bof');
            data = flipud(rot90(reshape(fread(fid,obj.Ny*obj.Nz, 'uint32',0,'l'),obj.Ny,obj.Nz)));
            fclose(fid);
                    
            % read data
            try
	            fid = fopen(filepath);
               fseek(fid, 3*512,'bof');
               data = flipud(rot90(reshape(fread(fid,obj.Ny*obj.Nz, 'uint32',0,'l'),obj.Ny,obj.Nz)));
               fclose(fid); 
            catch e
                warning(e.message);
                fprintf(1,'Debug information: \n');
                fprintf(1,'File path: %s\n',filepath);
                rethrow(e)
            end
        end
    end
end