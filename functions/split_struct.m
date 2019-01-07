function [varargout] = split_struct(struct_name, arglist)
% SPLIT_STRUCT  Description.
%   
%   result = split_struct(structure, cell_array_of_arguments) 
%
% The following arguments are supported:
%     structure:: [] (required)
%       A structure containing multiple fields.
%
%     cell_array_of_arguments:: [] (required)
%       A cell array of field names. If a field name was found in the
%       structure, it is added to the output argument list. 
%
% Example:
%   pca_result = pca(a_2d_diffraction_pattern,qy,qz,mask,[]);
%   [w,angle] = split_struct(pca_result,{'w','angle'});
%
% Output arguments:
%   This function outputs a variable number of arguments, based on the
%   number of query field names. 
%
%
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

    jj = 1;
    for ii = 1:numel(arglist)
        arg = arglist{ii};
        [isInStruct,val] = isfield_nested(struct_name,arg);
        if isInStruct
            varargout{jj} = val;
            jj = jj + 1;
        else
            warning([arg ' is not an element of input structure.']);
        end
    end
end