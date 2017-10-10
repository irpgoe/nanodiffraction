function out = sort_scan(data,rows,cols,type)

    % data input type
    if iscell(data)
        out = cell(rows,cols,size(data,2));
    elseif isstruct(data)
        out = repmat(data(1),rows,cols,size(data,2));
    else
        out = zeros(rows,cols,size(data,2));
    end

    % transform row to column
    if size(data,1) == 1
        data = data';
    end
    
    % treat columns seperately
    for ii = 1:size(data,2)
        data_per_col = data(:,ii);
        switch type
            case 'horz'
                data_per_col = transpose(reshape(data_per_col,[cols rows]));
            case 'horzalt'
                data_per_col = transpose(reshape(data_per_col,[cols rows]));
                data_per_col(2:2:end,:) = fliplr(data_per_col(2:2:end,:));
            case 'vert'
                data_per_col = reshape(data_per_col,[rows cols]);
            case 'vertalt'
                data_per_col = reshape(data_per_col,[rows cols]);
                data_per_col(:,2:2:end) = flipud(data_per_col(:,2:2:end));
            otherwise
                error('Choose between <horz>,<horzalt>,<vert>,<vertalt>');
        end
        out(:,:,ii) = data_per_col;
    end
end