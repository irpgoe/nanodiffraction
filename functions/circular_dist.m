function [circMean, circVar] = circular_dist(phi,I)
    % CIRCULAR_DIST  Calculates the circular mean and variance of a 
    % distribution based on 
    % https://en.wikipedia.org/wiki/Directional_statistics
    %
    % Accepted arguments:
    % phi, I, where I is usually given in counts and phi usually 
    % within the range [-180, 180] (because of central symmetry)

    % correct phi
    phi = phi*2-180;
    
    % number of data points
    N = numel(phi);
    
    % transform to complex number
    z = zeros(N,1);
    for jj = 1:N
        z(jj) = I(jj)*exp(1i*phi(jj)*pi/180);
    end

    % circular mean
    circMean = angle(sum(z)/N)*180/pi;
    
    % correct circ. mean
    circMean = (circMean+180)/2;
    
    % circular variance
    c = sum(cos(angle(z)));
    s = sum(sin(angle(z)));
    r = sqrt(c^2 + s^2)/N;
    circVar = 1-r; 
end    