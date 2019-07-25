function [Y, f] = Fourier(x, fs)

Fs = fs;                          % Sampling frequency
T = 1/Fs;                         % Sample time
L = length(x);                    % Length of signal
y=x;

NFFT = length(x); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
Y=Y(1:NFFT/2+1);
Y=2*(abs(Y).^2);
f = Fs/2* linspace(0,1,NFFT/2+1);

%Plot single-sided amplitude spectrum.
% figure
% semilogx(f,abs(Y),'color',[0 0 0]) 
% 
% xlabel('Frequency [Hz]')
% ylabel('FT Power')

end