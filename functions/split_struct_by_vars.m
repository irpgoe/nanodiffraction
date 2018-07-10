function [varargout] = split_struct_by_vars(struct_name, arglist)
    
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