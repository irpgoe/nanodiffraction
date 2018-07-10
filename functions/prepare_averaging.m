function [x,sortIndex,n,bin] = prepare_averaging(grid,mask,sel,bins)
    if isempty(bins)
        bins = 360;
    end

    % sort grid
    x = grid(~(mask | ~sel)); 
    [x, sortIndex] = sort(x);
    
    % binning of grid
    binEdge = linspace(min(x),max(x),bins);
    [n,~,bin] = histcounts(x,binEdge);
    
    % remove trailing zeros
    x = binEdge(1:sum(n~=0));
    
end