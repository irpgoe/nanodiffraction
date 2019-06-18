function [varargout] = split_struct(struct_name, arglist)
% SPLIT_STRUCT  Returns each indicated field as a separate variable.
%   
%   ``result = split_struct(structure, cell_array_of_arguments)``
%
% 
%
% Parameters
% ----------
%   structure: Structure
%       A structure that contains any number of fields. It can also be a nested structure
%
%   cell_array_of_arguments: Cell array
%       A cell array of field names. Each field name will be searched in the (possibly nested) structure. 
%       If a field name was found in the structure, it is added to the output argument list. 
%
% Returns
% -------
%   Variable-length output argument list. The length of the list is given by the number of fields that matched a query field name.
%
% Notes
% -----
% Example 1:
%
% .. code-block:: matlab
%
%       test_structure = struct('w',magic(25),'angle',pi/6);
%       [w,angle] = split_struct(test_structure,{'w','angle'});
%
% Example 2:
%
% .. code-block:: matlab
%
%       test_structure = struct('w', magic(25), 'angle', pi/6, struct('power', 2^3, 'prime', 13));
%       [w,angle,prime] = split_struct(test_structure,{'w','angle','prime'});
%

% Copyright 2017 Institute for X-ray Physics (University of Göttingen)
%
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
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