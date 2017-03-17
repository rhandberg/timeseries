function [Power,f] = fftps(t, x, full)
% Calculates the powerspectrum
% Input:
%    t - Time-vector in seconds
%    x - Data-vector

	% Number of data-points:
	N = length(t);
	T = t(end)-t(1); % max værdi - sluttid
	dT = T / N;

	% Pad zeros to next power of 2 to make FFT faster:
	%x = [x; zeros(1, 2.^nextpow2(N)-N)];
	%N = 2.^nextpow2(N);
	
	% Beregn Powerspektrum:
	Power = fft(x, N);
	Power = Power .* conj(Power); % Det samme som at tage længden i anden

	% Normér poweraksen:
	Power = Power ./ ((N^2)/4);
    
	% Lav frekvens-akse
	f = 0:(1/(dT*(N-1))):(1/dT);

	% Tag halvdelen:
    if nargin < 3 || ~strcmp(full,'full')
        f = f(1:ceil(N/2));
        Power = Power(1:ceil(N/2))';
    end;