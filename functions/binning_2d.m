function image = binning_2d(image,binx,biny)                
    % BINNING_2D  takes in an image and bins it according to the binning
    % factors BINX and BINY.
    %
    %   Usage:
    %       BINNED_IMAGE = BINNING_2D(IMAGE, BINX, BINY)
    %
    %   Arguments:
    %       image:: [no default]
    %           2d array.
    %
    %       binx:: [no default]
    %           Integer binning factor along the horizontal dimension.
    %
    %       biny:: [no default]
    %           Integer binning factor along the vertical dimension.

    [m,n] = size(image);
    nnew = n - mod(n,binx);
    mnew = m - mod(m,biny);

    % now it is divisible by binx, biny
    image = image(1:mnew,1:nnew);

    % binning along z (correct)
    image=sum( reshape(image,biny,[]),1); 
    image=reshape(image,mnew/biny,[]); 

    % binning along y
    image=sum(reshape(image',binx,[]),1); 
    image=reshape(image,nnew/binx,[])';
    
end