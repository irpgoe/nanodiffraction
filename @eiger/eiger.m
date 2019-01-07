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
        c;
        fpf;
        h5config;
    end
    methods
        
        function obj = eiger(config,fpf)
            obj.fpf = fpf;
            obj.c = config;
            obj.c.fast_comp = true;
            
            % default h5 configuration
            obj.h5config.master_file = config.path;
            obj.h5config.data_entry = 'data_000001';
            obj.h5config.initialized = 0;
            
            obj.repair_frames_per_file();
            
            fprintf(1,'Initialized Eiger module with the following settings:\n');
            fprintf(1,'Data path: %s\n',obj.c.path);
            fprintf(1,'Frames per file: %d\n',obj.fpf);
        end
        
        function data = read_standard(obj,fn)
            index = floor(double(fn-1)/double(obj.fpf))+1;
            remainder = mod(double(fn-1),double(obj.fpf))+1;
            data_entry = ['/entry/data/data_' sprintf('%06i',index)];
            file_in_data_entry = double(remainder);
            % read data
            try
                data = obj.eiger_read(obj.c.path,data_entry,[1 1 file_in_data_entry],[2070 2167 1]); 
%                 data = flipud(permute(flip(data,2),[2 1 []]));                 
                data = transpose(data);
            catch e
                warning(e.message);
                fprintf(1,'Debug information: \n');
                fprintf(1,'Static path: %s\nData Entry:%s\nFile in data entry:%d\nFrames per file:%d\n',obj.c.path,data_entry,file_in_data_entry,obj.fpf);
                rethrow(e)
            end
        end
                
        function data = read(obj,fn)
            index = floor(double(fn-1)/double(obj.fpf))+1;
            remainder = mod(double(fn-1),double(obj.fpf))+1;
            data_entry = sprintf('data_%06i',index);
            file_in_data_entry = double(remainder);
            
            % check for updates
%             [upd_master, upd_data_entry] = update_h5config(obj,obj.c.path,data_entry);
            
            % read data
            try
                % Normally, I would call h5read() here, but now, lets use
                % low-level functions instead
                
%                 if obj.c.fast_comp && ~(upd_master || upd_data_entry) && obj.h5config.initialized
% 
%                     [file_space_id,dset_id,type_id,mem_space_id,block] = split_struct(obj.h5config,{'file_space_id','dset_id','type_id','mem_space_id','block'});
%                     offset = [file_in_data_entry 0 0];
%                     H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',offset,[],[],block);
% 
%                     % read data
%                     data = H5D.read(dset_id,type_id,mem_space_id,file_space_id,'H5P_DEFAULT');
%                 else 
                    
                    % open an existing master file
                    fid = H5F.open(obj.c.path,'H5F_ACC_RDONLY','H5P_DEFAULT'); % obj.c.path = /.../myfile_master.h5
                    
                    % open an existing group
                    gid = H5G.open(fid,'/entry/data','H5P_DEFAULT');
                    
                    % open an existing dataset
                    dset_id = H5D.open(gid,data_entry);

                    % what is the required space?
                    file_space_id = H5D.get_space(dset_id);
                    
                    % Gets the dimensions of the data set including its rank. 
                    % The HDF5 library uses C-style ordering for 
                    % multidimensional arrays, while MATLAB® uses FORTRAN-style
                    % ordering. The h5_dims and h5_maxdims assume C-style 
                    % ordering.
                    [ndims,h5_dims] = H5S.get_simple_extent_dims(file_space_id);
                    
                    % what is the type of the data set?
                    type_id = H5D.get_type(dset_id);

                    % create memory space for a single data frame
                    block = [1 h5_dims(2) h5_dims(3)]; 
                    mem_space_id = H5S.create_simple(ndims,block,[]);

                    % selecting an hyperslab (single frame)
                    offset = [file_in_data_entry-1 0 0]; 
                    H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',offset,[],[],block);

                    % read data
                    data = H5D.read(dset_id,type_id,mem_space_id,file_space_id,'H5P_DEFAULT');

                    % close everything
                    H5T.close(type_id);
                    H5D.close(dset_id);
                    H5F.close(fid);
%                 end
                
                % uint16 to double and row-first to column-first indexing
%                 data = double(data);
                data = transpose(double(data));
            catch e
                warning(e.message);
                fprintf(1,'Debug information: \n');
                fprintf(1,'Static path: %s\nData Entry:%s\nFile in data entry:%d\nFrames per file:%d\n',obj.c.path,data_entry,file_in_data_entry,obj.fpf);
                rethrow(e)
            end
        end
        
        function init_h5config(obj)
            master_file = obj.c.path;
            data_entry = 'data_000001';
            
            % open an existing master file
            fid = H5F.open(master_file,'H5F_ACC_RDONLY','H5P_DEFAULT'); % obj.c.path = /.../myfile_master.h5
            % open an existing group
            gid = H5G.open(fid,'/entry/data','H5P_DEFAULT'); 
            % open an existing dataset
            dset_id = H5D.open(gid,data_entry);

            % what is the required space?
            file_space_id = H5D.get_space(dset_id);
            % Gets the dimensions of the data set including its rank. 
            % The HDF5 library uses C-style ordering for 
            % multidimensional arrays, while MATLAB® uses FORTRAN-style
            % ordering. The h5_dims and h5_maxdims assume C-style 
            % ordering.
            [ndims,h5_dims] = H5S.get_simple_extent_dims(file_space_id);
            % In the case of an Eiger 4M, this should be [2070 2167
            % 2000] where the last number can be arbitrary
            % matlab_dims = fliplr(h5_dims);
            % what is the type of the data set?
            type_id = H5D.get_type(dset_id);

            % create memory space for a single data frame
            block = [1 h5_dims(2) h5_dims(3)];
            mem_space_id = H5S.create_simple(ndims,block,[]);

            % save configuration
            obj.h5config = struct('fid',fid,'gid',gid,'dset_id',dset_id,...
                'file_space_id',file_space_id,'ndims',ndims,'h5_dims',h5_dims,...
                'type_id',type_id,'block',block,'mem_space_id',mem_space_id,...
                'initialized',1,'master_file',master_file,'data_entry',data_entry);
        end
        
        function close_eiger_files(obj)
            % close everything
            H5T.close(obj.h5config.type_id);
            H5D.close(obj.h5config.dset_id);
            H5F.close(obj.h5config.fid);
        end
        
        function [upd_master, upd_data_entry] = update_h5config(obj,master_file,data_entry)
            
            % default
            upd_master = 1;
            upd_data_entry = 1;
            
            % did the master file change?
            if strcmp(master_file,obj.h5config.master_file)
                % no change
                upd_master = 0;
            end
            
            % did the data entry change?
            if strcmp(data_entry,obj.h5config.data_entry)
                % no change
                upd_data_entry = 0;
            end                        
        end
        
        function repair_frames_per_file(obj)
            % REPAIR_FRAMES_PER_FILE  calculates frames per file when using 
            % the h5 container format
            %
            %   REPAIR_FRAMES_PER_FILE()
            
            % starting point
            n = 1;
            m = 100;
            fpf = 1;
            cond = 1;
            while cond
                
                % try to read
                try 
                    data_entry = ['/entry/data/data_' sprintf('%06i',1)];
                    obj.eiger_read(obj.c.path,data_entry,[1 1 n],[2070 2167 1]);
                catch
                    % we have gone out of bounds
                    n = n - m; % go back to value before
                    m = m / 10; % reduce sampling rate
                    if m < 1 
                        cond = 0;
                    end
                end
                % last working condition
                fpf = n;
                % we are still within bounds, hence step forward
                n = n + m;
                
            end
            obj.fpf = fpf;
        end
    end
end