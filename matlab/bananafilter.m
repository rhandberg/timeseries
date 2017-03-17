function [t2, D] = bananafilter(t, x, B1, factor)
%Perform the Banana-filter on timeseres.
%
%Input:
%   t      - Vector of times.
%   x      - Vector of data.
%   B1     - Width of central part of filter.
%   factor - Factor between width of wings, B2, and B1.
%            Default = 2
%Returns:
%   t2     - Vector of modified times.
%   D      - Vector of filtered data.
%

    if (nargin < 4)
        factor = 2;
    end;

    N = length(t);
    tmin = min(t);
    dt = (max(t)-tmin)/N;
    B2 = factor*B1;
    D = zeros(N, 1);
    t2 = zeros(N, 1);

    for i = 1:N
        ti = (i-1)*dt+tmin;
        
        mask = find(t > ti-B1/2-B2 & t < ti-B1/2);
        if isempty(mask)
            W1 = 0;
        else
            W1 = mean( x(mask) );
        end;
        
        mask = find(t < ti+B1/2+B2 & t > ti+B1/2);
        if isempty(mask)
            W2 = 0;
        else
            W2 = mean( x(mask) );
        end;

        C  = mean( x(t >= ti-B1/2 & t <= ti+B1/2) );
        
        D(i) = C - (W1 + W2)/2;
        t2(i) = ti;
    end;
