function [isFieldResult,val] = isfield_nested(inStruct, fieldName)
% ISFIELD_NESTED checks wether a field exists within a structure
%
%   ``[isField,value] = isfield_nested(inStruct, fieldName)
%
% Parameters
% ----------
%   inStruct: Structure or structure array
%       The name of the structure or an array of structures to search
%
%   fieldName: String 
%       The name of the field for which the function searches
%
% Returns
% -------
%   isField: boolean
%       Boolean value; true = field exists, false = field does not exist
%
%   value: numeric
%       Value of the identified field. If field was not found inside the structure, a value of NaN is returned.
%
% Notes
% -----
% This function was obtained from: 
% https://de.mathworks.com/matlabcentral/answers/95923-is-there-a-matlab-function-that-can-check-if-a-field-exists-in-a-matlab-structure

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

    isFieldResult = 0;
    f = fieldnames(inStruct(1));
    val = NaN;
    for i=1:length(f)
        if(strcmp(f{i},strtrim(fieldName)))
            isFieldResult = 1;
            val = inStruct(1).(f{i});
            return;
        elseif isstruct(inStruct(1).(f{i}))
            [isFieldResult,val] = isfield_nested(inStruct(1).(f{i}), fieldName);
            if isFieldResult
                return;
        end
    end
end