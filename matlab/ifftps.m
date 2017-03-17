function z = ifftps(P)
% Beregner signalet i tidsdomænet ud fra powerspektret.

    % Antallet af datapunkter:
    N = length(P);

    % Brug invers FFT til at omregne til tids-domænet:
    z = ifft(sqrt(P), 'symmetric'); % sqrt fordi ifft åbenbart skal have amp-spektrum som input
    z = real(z); % For at fjerne afrundingsfejl
    z = z * N/2; % Mere syg normering