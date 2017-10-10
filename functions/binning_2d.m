function image = binning_2d(image,binx,biny)                

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