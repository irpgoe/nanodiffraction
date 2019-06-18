function config = update_defaults(defaults,opts)
% UPDATE_DEFAULTS  overwrites values in a structure with given values from
% a reference structure
%
%   ``config = update_defaults(defaults,opts)``
%
% Parameters
% ----------
% defaults: Structure
%     Structure containing default entries.
%
% opts: Structure
%     Structure containing entries that should be used to overwrite
%     entries in the default structure.
%
% Returns
% -------
% config: Structure
%   A structure that contains all fields from the default structure
%
% Notes
% -----
% Example 1:
%
% .. code-block:: matlab
%
%   config = update_defaults(struct('a',1,'b','test'),struct('b','replaced'));
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

    
    fields = fieldnames(opts);
    for f = 1:numel(fields)
        if ~isfield(defaults,fields{f})
            warning(sprintf('You provided an illegal option: %s',fields{f}));
        else
            defaults.(fields{f}) = opts.(fields{f});
        end
    end
    config = defaults;

end