function data = baseline(data, winSize, iterations)
    % BASELINE  performs...
    %   
    %   [RESULT] = BASELINE(DATA,WINSIZE,ITERATIONS) 
    %
    %   The following options are required:
    %
    %     DATA::
    %       The data that will be processed.
    %
    %     WINSIZE::
    %       Integer > 0
    %
    %     ITERATIONS::
    %       Integer > 0
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
    
    
    m = data;
    
    for ii = (1+winSize):(numel(data)-winSize)
        m(ii) = 0.5*(data(ii-winSize) + data(ii+winSize));
    end
    sel = m < data;
    data(sel) = m(sel);

    if iterations == 0
        return;
    end
    data = baseline(data, winSize, iterations-1);
end