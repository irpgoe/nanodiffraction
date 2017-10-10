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
classdef eiger<handle
    properties
        eiger_read = @(str1,str2,vec1,vec2) double(h5read( str1, str2, vec1, vec2 ));
        static_path;
        framesPerFile;
    end
    methods
        function obj = eiger(beamline,prepath,newfile,eigerNr,framesPerFile,slash)
            obj.framesPerFile = framesPerFile;
            switch beamline
                case 'p10'
                    obj.static_path = [slash prepath...
                        slash 'detectors'...
                        slash 'eiger'...
                        slash newfile...
                        slash newfile '_'...
                        sprintf('%05i',eigerNr) '_master.h5'];
                case 'id13'
                    obj.static_path = [slash prepath...
                        slash newfile '_'...
                        num2str(eigerNr) '_master.h5'];
                case 'none'
                    obj.static_path = [slash prepath...
                        slash newfile '_'...
                        sprintf('%05i',eigerNr) '_master.h5'];
            end
        end
        function data = read(obj,fn)
            index = floor(double(fn-1)/double(obj.framesPerFile))+1;
            remainder = mod(double(fn-1),double(obj.framesPerFile))+1;
            data_entry = ['/entry/data/data_' sprintf('%06i',index)];
            file_in_data_entry = double(remainder);
            % read data
            try
                data = flipud(permute(flip(obj.eiger_read(obj.static_path,data_entry,[1 1 file_in_data_entry],[2070 2167 1]),2),[2 1 []])); 
            catch e
                warning(e.message);
                fprintf(1,'Debug information: \n');
                fprintf(1,'Static path: %s\nData Entry:%s\nFile in data entry:%d\nFrames per file:%d\n',obj.static_path,data_entry,file_in_data_entry,obj.framesPerFile);
                rethrow(e)
            end
        end
    end
end