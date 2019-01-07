function save_progress(script_base,script_fn,location,desc)
% SAVE_PROGRESS  saves ... 
%   
%   save_progress(script_base,script_fn,location,desc) 
%
% The following arguments are supported:
%     script_base:: ['./experiments'] (required)
%       The data that will be processed. If script_base is empty, then the
%       default value will be used.
%
%     script_fn:: [] (required)
%       Name of the script that is currently being used, e.g.ls2522
%
%     location:: [] (required)
%       location where result will be saved. Note, that the location should
%       not contain an ending slash.
%
%     desc:: [] (required)
%       String that describes what you did.
%
% Example:
%   save_progress([],ls2728,'/results/ls2728','I added a few things');
%
% Output arguments:
%   This function does not return any arguments.
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

    if isempty(script_base)
        script_base = './experiments';
    end
    
    timestamp = datestr(now,'yyyy-mm-dd_HH:MM');
    syscall = ['cp ' script_base '/' script_fn '.m ' location '/' timestamp '.m'];
    expression = ['save(''', [location '/' timestamp '.mat'], ''',','''-v7.3'')'];

    fprintf(1,'\n');
    disp('You are currently in the following directory:')
    disp(pwd)
    fprintf(1,'\n');
    disp('Workspace saving:')
    disp(expression);
    fprintf(1,'\n');
    disp('Script backup:')
    disp(syscall);
    fprintf(1,'\n');
    disp('Appending logbook: ');
    disp(['./logs/' script_fn '.dat']);
    fprintf(1,'\n');
    answer = input('Do you wish to proceed? (y|n) ','s');
    
    if strcmp(answer,'y')
        disp('executing...')
    
        fid = fopen(['./logs/' script_fn '.dat'],'a');
        fprintf(fid,'Current time: %s\n%s\n\n',timestamp,desc);
        fclose(fid);
        
        evalin('base',expression);

        system(syscall);
    
    else
        disp('Nothing done.')
    end
    
    % append to description

end
