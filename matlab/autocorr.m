function ac = autocorr(x, kspan)
%Calculate the Autocorrelation of x.
%
% Input:
%     x     - Vector of lenght N for which to calculate the autocorrelation.
%     kspan - Optional. Span in k to calculate the autocorr.
%             Default = 1:N-1
% Returns:
%    ac  - Vector of length N-1 with the autocorrelation of x.
%
% Example: Finding the large seperation in spectrum.
%    A = spec(timeseries);
%    ac = autocorr(A);
%    nus = nu(2:end);
%    plot(nus, ac);
%
% See also SPEC.

	N = length(x);
	if nargin < 2
		kspan = 1:N-1;
	end;
    ac = zeros(size(kspan));

    for k = kspan
        % Calculate the mean-values:
        xmean = mean(x(1:N-k));
        ymean = mean(x(1+k:N));

        % Calculate the autocorrelation:
        i = 1:N-k;
        s = sum( (x(i)-xmean) .* (x(i+k)-ymean) );
        
        ac(k) = s/(N-k);
    end;