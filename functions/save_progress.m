function save_progress(script_fn,location,desc)
% SAVE_PROGRESS  saves ... 
%   Syntax: SAVE_PROGRESS(SCRIPT_FN,LOCATION,DESC)
%   script_fn:: Name of the script that is currently being used, e.g.ls2522
%   location:: location where result will be saved
%   desc:: String that describes what you did
%   
%   Note, that the location should not contain an ending slash
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

    timestamp = datestr(now,'yyyy-mm-dd_HH:MM');
    syscall = ['cp ./experiments/' script_fn '.m ' location '/' timestamp '.m'];
    expression = ['save(''', [location '/' timestamp '.mat'], ''',','''-v7.3'')'];

    disp('Workspace saving:')
    disp(expression);
    disp('')
    disp('Script backup:')
    disp(syscall);
    disp('')
    disp('Appending logbook: ');
    disp([script_fn '.dat']);
    disp('')    
    answer = input('Do you wish to proceed? (y|n)','s');
    
    if strcmp(answer,'y')
        disp('executing...')
    
        fid = fopen([script_fn '.dat'],'a');
        fprintf(fid,'Current time: %s\n%s\n\n',timestamp,desc);
        fclose(fid);
        
        evalin('base',expression);

        system(syscall);
    
    else
        disp('Nothing done.')
    end
    
    % append to description

end