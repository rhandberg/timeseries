function z = ifftps(P)
% Beregner signalet i tidsdom�net ud fra powerspektret.

    % Antallet af datapunkter:
    N = length(P);

    % Brug invers FFT til at omregne til tids-dom�net:
    z = ifft(sqrt(P), 'symmetric'); % sqrt fordi ifft �benbart skal have amp-spektrum som input
    z = real(z); % For at fjerne afrundingsfejl
    z = z * N/2; % Mere syg normering