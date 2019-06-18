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
classdef xia<handle
    properties
        fpf;
        path;
        filename;
    end
    methods
        function obj = xia(config)
                        
            defaults = struct('path','',...
                'fpf',128,...
                'filename',@(pre,n) [pre sprintf('%04i',n) '.edf']);
            config = update_defaults(defaults,config);

            % get configuration
            [obj.path, obj.fpf, obj.filename] = split_struct(config,{'path','fpf','filename'});
            
        end
        
        function data = read(obj,fn)

	        % file numbers start with 0
	        remainder = mod(double(fn-1),double(obj.fpf))+1;
	        filename = obj.filename(obj.path, floor((fn-1)/obj.fpf));

            % read data
            try
                tmp = read_spectrum(filename);
	           	data = tmp.data(:,double(remainder)); 
            catch e
                warning(e.message);
                fprintf(1,'Debug information: \n');
                fprintf(1,'File path: %s\n',filename);
                rethrow(e)
            end
        end
        
        function write(obj,data)
            
        end
    end
end