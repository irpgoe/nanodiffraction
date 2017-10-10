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
classdef rayonix<handle
    properties
        prepath;
        suffix;
    end
    methods
        function obj = rayonix(beamline,prepath,newfile,slash)
				obj.suffix = '_raw.edf';            
            switch beamline
                case 'id02'
                    obj.prepath = [slash prepath...
                        slash newfile '_'];
            end
        end
        function data = read(obj,fn)
        
            filepath = [obj.prepath sprintf('%05i',fn) obj.suffix];	
%                     fprintf(1,'argument 1: %s\n',filepath{1});
            fid = fopen(filepath);
            fseek(fid, 8*512,'bof');
            data = flipud(rot90(reshape(fread(fid,960*960, 'uint16',0,'l'),960,960)));
            fclose(fid);
                    
                    
            % read data
            try
	            fid = fopen(filepath);
               fseek(fid, 8*512,'bof');
               data = flipud(rot90(reshape(fread(fid,960*960, 'uint16',0,'l'),960,960)));
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