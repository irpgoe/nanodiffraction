function arr = cellmat2array(cellmat)

    [sx,sy] = size(cellmat{1});
    arr = reshape(cell2mat(cellmat),sx,sy,[]);
end