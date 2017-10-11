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

classdef pilatus<handle
    properties
        static_path;
        suffix;
    end
    methods
        function obj = pilatus(beamline,prepath,newfile,slash)
            
            obj.suffix = '.cbf';            
            if ~ispc
                prepath = [slash prepath];
            end
            
            switch beamline
                case 'p10'
                    obj.static_path = [prepath...
                        slash 'detectors'...
                        slash 'p300'...
                        slash newfile...
                        slash newfile '_'];
            end
        end
        function data = read(obj,fn)
        
            filepath = [obj.static_path sprintf('%05i',fn) obj.suffix];
        		
            % read data
            try
                data = flipud(read_pilatus(filepath,'p300')); 
            catch e
                warning(e.message);
                fprintf(1,'Debug information: \n');
                fprintf(1,'File path: %s\n',filepath);
                rethrow(e)
            end
        end
    end
end