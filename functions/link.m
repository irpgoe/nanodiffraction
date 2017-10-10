function link(classA,classB,varargin)
% LINK  links class A to class B
%
%   LINK(CLASSA,CLASSB,VARARGIN)
%
%   The following arguments are required:
%       CLASSA:: ()
%           Class that will be linked from.
%
%       CLASSA:: ()
%           Class that will be linked to.
%
%       VARARGIN:: ()
%           Additional class pairs can be given.
%
% These combinations are possible
%   fpath > nanodiffraction (tells nanodiffraction where the data is)
%   nanodiffraction > visualization (tells visualization how the
%   experiment was done)
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
    
    if isa(classA,'files')
        % fpath class is linkable to experiment and display
        if isa(classB,'nanodiffraction')
            classB.data = classA;
            disp('Linked data to experiment');
        elseif isa(classB,'display')
        end
    elseif isa(classA,'nanodiffraction')
        % fpath class is linkable to experiment and display
        if isa(classB,'nanodiffraction')
        elseif isa(classB,'display')
            classB.exp = classA;
            disp('Linked experiment to visualization');
        end
    elseif isa(classA,'display')
        % fpath class is linkable to experiment and display
        if isa(classB,'nanodiffraction')
            classB.vis = classA;
        elseif isa(classB,'display')
        end
    else
        warning('You seem to be trying to link classes together, that are unlinkable. Please type "help link" for usage');
    end
    
    if nargin == 3
        warning('Uneven number of classes. Please type "help link" for usage');
    end
    
    if nargin > 3
        % if more classes are to be linked, then call link recursively
        link(varargin{1},varargin{2},varargin{3:end})
    end
end