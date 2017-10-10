function [frame,vararg_remain] = read_spectrum(filename,varargin)
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

% 0: no debug information
% 1: some feedback
% 2: a lot of information
debug_level = 0;

% initialize return argument
frame = struct('header',[], 'data',[]);

% check minimum number of input arguments
if (nargin < 1)
    image_read_sub_help(mfilename,'edf');
    error('At least the filename has to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 2)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 1 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 1)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end

    
% set default values for the variable input arguments and parse the named
% parameters: 
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        otherwise
            % pass further arguments on to fopen_until_exists
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end


% expected maximum length for the text header
max_header_length = 4096;


% try to open the data file
if (debug_level >= 1)
    fprintf('Opening %s.\n',filename);
end
[fid,vararg_remain] = fopen_until_exists(filename,vararg);
if (fid < 0)
    return;
end

% read all data at once
[fdat,fcount] = fread(fid,'uint8=>uint8');

% close input data file
fclose(fid);
if (debug_level >= 2)
    fprintf('%d data bytes read\n',fcount);
end

% search for end of header signature within the expected maximum length of
% a header
end_of_header_pos = 512;
max_pos = min( max_header_length,length(fdat) );

while ((end_of_header_pos < max_pos) && (fdat(end_of_header_pos-1) ~= '}'))
    end_of_header_pos = end_of_header_pos +512;
end
if (end_of_header_pos >= max_pos)
    disp('no header end signature found');
end
if (debug_level >= 2)
    fprintf('Header length is %d bytes.\n',end_of_header_pos);
end
data_length = fcount - end_of_header_pos;

% convert the header to lines of a cell array
frame.header = char_to_cellstr( char(fdat(1:(end_of_header_pos-1))') );

% check for opening parenthesis
if (frame.header{1} ~= '{')
    error([filename ': EDF start ,''{'' not found in first line ''' ...
        frame.header{1} '''' ]);
end


% extract the mandatory information for data extraction from the header:
byte_order = get_hdr_val(frame.header,'ByteOrder',' = %s',1);
dim1 = get_hdr_val(frame.header,'Dim_1',' = %d',1);
dim2 = get_hdr_val(frame.header,'Dim_2',' = %d',1);
data_type = get_hdr_val(frame.header,'DataType',' = %s',1);
if (debug_level >= 2)
    fprintf('Byte order is %s\n',byte_order);
    fprintf('Frame dimensions are %d x %d.\n',dim2,dim1);
    fprintf('Data type is %s\n',data_type);
end

% determine number of bytes per pixel
switch data_type
    case 'UnsignedByte',
        bytes_per_pixel = 1;
        data_class = 'uint8';
    case 'UnsignedShort',
        bytes_per_pixel = 2;
        data_class = 'uint16';
    case {'SignedInteger','UnsignedInteger','UnsignedInt','UnsignedLong'}
        bytes_per_pixel = 4;
        data_class = 'uint32';
    case {'Float','FloatValue','Real'}
        bytes_per_pixel = 4;
        data_class = 'single';
    case 'DoubleValue'
        bytes_per_pixel = 8;
        data_class = 'double';
    otherwise
        error('unsupported data type %s',data_type);
end
no_of_bytes = bytes_per_pixel * dim1 * dim2;
if (debug_level >= 2)
    fprintf('%d bytes per pixel, %d in total expected, %d available\n',...
        bytes_per_pixel,no_of_bytes,data_length);
end

% check length of available data
if (no_of_bytes > data_length)
    error('%d data bytes expected, %d are available',...
        no_of_bytes,data_length);
end

% compare file with machine byte order, swap if necessary
[str,maxsize,endian] = computer;
if (((strcmp(byte_order,'HighByteFirst')) && (endian == 'L')) || ...
    ((strcmp(byte_order,'LowByteFirst')) && (endian == 'H')))
    if (debug_level >= 2)
        fprintf('Machine byte order is %s: swapping data bytes\n',...
            endian,bytes_per_pixel);
    end
    dat = fdat(end_of_header_pos+1:end_of_header_pos+no_of_bytes);
    dat = reshape(dat,bytes_per_pixel,[]);
    dat = flipud(dat);
    fdat(end_of_header_pos+1:end_of_header_pos+no_of_bytes) = dat(:);
end

% extract the frame from the binary data
[frame.data] = ...
    double(reshape(typecast(fdat(end_of_header_pos+1:end_of_header_pos+no_of_bytes),...
                            data_class),...
                   dim1,dim2));

end




% Filename: $RCSfile: get_hdr_val.m,v $
%
% $Revision: 1.3 $  $Date: 2008/08/28 18:47:31 $
% $Author: bunk $
% $Tag: $
%
% Description:
% Find text signature in a bunch of cell strings from a file header and
% return the following value in the specified format. Example:
% no_of_bin_bytes = get_hdr_val(header,'X-Binary-Size:','%f',1);
% The last parameter specifies whether the macro should exit with an error
% message if the text signature has not been found.
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies:
% none
%
%
% history:
%
% May 7th 2008, Oliver Bunk:
% add number of input argument check and brief help text
%
% April 25th 2008, Oliver Bunk: 1st version
%
function [outval,line_number,err] = get_hdr_val(header,signature,format,...
    exit_if_not_found)

% initialize output arguments
outval = 0;
line_number = 0;
err = 0;

if (nargin ~= 4)
    fprintf('Usage:\n');
    fprintf('[value,line_number,error]=%s(header,signature,format,exit_if_not_found);\n',...
        mfilename);
    fprintf('header             cell array with text lines as returned by cbfread or ebfread\n');
    fprintf('signature          string to be searched for in the header\n');
    fprintf('format             printf-like format specifier for the interpretation of the value that follows the signature\n');
    fprintf('exit_if_not_found  exit with an error in case either the signature or the value have not been found\n');
    error('Wrong number of input arguments.\n');
end

% search for the signature string
pos_found = strfind(header,signature);

% for sscanf the percentage sign has a special meaning
signature_sscanf = strrep(signature,'%','%%');

% loop over the search results for all header lines
for (ind=1:length(pos_found))
    % if the signature string has been found in this line
    if (length(pos_found{ind}) > 0)
        % get the following value in the specified format
        [outval,count] = sscanf(header{ind}(pos_found{ind}:end),...
            [signature_sscanf format]);
        % return an error if the signature and value combination has not
        % been found (i.e., the format specification did not match)
        if (count < 1)
            outval = 0;
            err = 1;
        else
            % return the first occurrence if more than one has been found
            if (count > 1)
                outval = outval(1);
            end
            % return the line number
            line_number = ind;
            return;
        end
    end
end

% no occurrence found
err = 1;
if (exit_if_not_found)
    error(['no header line with signature ''' signature ''' and format ' ...
        format ' found']);
end

return;
end




% May 9th 2008, Oliver Bunk: 1st version
%
function [fid,vararg_remain] = fopen_until_exists(filename,varargin)

if (nargin < 1)
    fprintf('Usage:\n');
    fprintf('[fid] = %s(filename [[,<name>,<value>],...]);\n',...
        mfilename);
    fprintf('filename                             name of the file to open\n');
    fprintf('The optional name value pairs are:\n');
    fprintf('''RetryReadSleep'',<seconds>           if greater than zero retry opening after this time (default: 0.0)\n');
    fprintf('''RetryReadMax'',<0-...>               maximum no. of retries, 0 for infinity (default: 0)\n');
    fprintf('''MessageIfNotFound'',<0-no,1-yes>     display a mesage if not found, 1-yes is default\n');
    fprintf('''ErrorIfNotFound'',<0-no,1-yes>       exit with an error if not found, default is 1-yes\n');
    fprintf('The file ID of the opened file is returned or -1 in case of failure.\n');
    error('Invalid number of input parameters.');
end

% check minimum number of input arguments
if (nargin < 1)
    display_help();
    error('At least the filename has to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 2)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 1 + length(varargin);
        end
    end
end
% check number of input arguments
if (rem(no_of_in_arg,2) ~= 1)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end

% set default values for the variable input arguments and parse the named
% parameters:

% If the file has not been found and if this value is greater than 0.0 than
% sleep for the specified time in seconds and retry reading the file.
% This is repeated until the file has been successfully read
% (retry_read_max=0) or until the maximum number of iterations is exceeded
% (retry_read_max>0).
retry_read_sleep_sec = 0.0;
retry_read_max = 0;

% exit with error message if the file has not been found
error_if_not_found = 1;

% display a message once in case opening failed
message_if_not_found = 1;

% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'RetryReadSleep'
            retry_read_sleep_sec = value;
        case 'RetryReadMax'
            retry_read_max = value;
        case 'MessageIfNotFound'
            message_if_not_found = value;
        case 'ErrorIfNotFound'
            error_if_not_found = value;
        otherwise
            vararg_remain{end+1} = name;
            vararg_remain{end+1} = value;
    end
end


% try to access the file entry
file_non_empty = 0;
dir_entry = dir(filename);

% if it has not been found or if it is empty
if ((isempty(dir_entry)) || (size(dir_entry,1) == 0) || ...
        (dir_entry.bytes <= 0))
    if (message_if_not_found)
        if (isempty(dir_entry))
            fprintf('%s not found',filename);
        else
            fprintf('%s found but of zero length',filename);
        end
    end
    % retry, if this has been specified
    if (retry_read_sleep_sec > 0.0)
        if (message_if_not_found)
            fprintf(', retrying\n');
        end
        % repeat until found or the specified number of repeats has been
        % exceeded (zero repeats means repeat endlessly)
        retry_read_ct = retry_read_max;
        while ((~file_non_empty) && ...
                ((retry_read_max <= 0) || (retry_read_ct > 0)))
            pause(retry_read_sleep_sec);
            dir_entry = dir(filename);
            if ((~isempty(dir_entry)) && (dir_entry.bytes > 0))
                file_non_empty = 1;
                % !!! arbitrary pause to increase the likelihood that an
                % HDF5 file is written !!!
                pause(1);
            end
            retry_read_ct = retry_read_ct -1;
        end
    else
        fprintf('\n');
    end
else
    file_non_empty = 1;
end

% open the file for read access
if (file_non_empty)
    fid = fopen(filename,'r');
else
    fid = -1;
end

% exit with an error message, if this has been specified and if the file
% could not be opened
if (fid < 0)
    if (error_if_not_found)
        error('file ''%s'' not found',filename);
    end
end
end








%
% Filename: $RCSfile: char_to_cellstr.m,v $
%
% $Revision: 1.3 $  $Date: 2008/06/23 13:49:50 $
% $Author: bunk $
% $Tag: $
%
% Description:
% Convert an array of text to a cell array of lines.
%
% Note:
% Used for making file headers accessible. 
%
% Dependencies:
% none
%
%
% history:
%
% June 22nd 2008, Oliver Bunk: bug fix for fliread adding the nl_only 
% parameter, to be replaced by named parameter later on
%
% May 9th 2008, Oliver Bunk: 1st version
%
function [outstr] = char_to_cellstr(inchars,nl_only)

if (nargin < 2)
    nl_only = 0;
end

% get positions of end-of-line signatures
eol_ind = regexp(inchars,'\r\n');
eol_offs = 1;
if ((length(eol_ind) < 1) || (nl_only))
    eol_ind = regexp(inchars,'\n');
    eol_offs = 0;
end
if (length(eol_ind) < 1)
    eol_ind = length(inchars) +1;
end
if (length(eol_ind) < 1)
    outstr = [];
    return;
end

% dimension return array with number of lines
outstr = cell(length(eol_ind),1);

% copy the lines to the return array, suppressing empty lines
start_pos = 1;
ind_out = 1;
for (ind = 1:length(eol_ind))
    end_pos = eol_ind(ind) -1;
    % cut off trailing spaces
    while ((end_pos >= start_pos) && (inchars(end_pos) == ' '))
        end_pos = end_pos -1;
    end
    % store non-empty strings
    if (end_pos >= start_pos)
        outstr{ind_out} = inchars(start_pos:end_pos);
        ind_out = ind_out +1;
    end
    
    start_pos = eol_ind(ind) +1 + eol_offs;
    ind = ind +1;
end

% resize cell array in case of empty lines
if (ind_out <= length(eol_ind))
    outstr = outstr(1:(ind_out-1));
end
end
