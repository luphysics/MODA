
function wt = wtAtfMod(sig, fs, fr, varargin)
% author:   Aleksandra Pidde, a.pidde@lancaster.ac.uk, aleksandra.pidde@gmail.com
% date:         30.12.2017
% update:       8.01.2018
%
% calculating wavelets at certain frequency
% INPUT
% sig:          signal
% fs:           sampling frequency [Hz]
% fr:           frequency for which wavelets are calculated
% varargin (optional): 
% 'f0':         tradeoff between the time and frequency resolutions, the higher, the
%               better frequency resolution, the lower the better time resolution,
%               default 1
% 'wavelet':    type of wavelet, {'morlet', 'modMorlet'}, default 'morlet'
% 'cutEdges':   cutting of the cone of influence, {true, false}, default
%               true
% 'd':          Gaussian decay in Morlet wavelet, according to
%               correction, default 2
% 'p':          normalization factor (Kaiser, A friendly guide to wavelets,
%               params from eq 3.6), default 0.5
%
% OUTPUT
% WT:           wavelet transform
% freq:         frequency vector

wavelet = 'morlet';
f0 = 1; 
cutEdges = true;
d = 2; 
p = 0.5;
if nargin >= 2 + 3
    for i = 1 : 2 : nargin - 3
        switch varargin{i}
            case 'wavelet'
                wavelet = varargin{i + 1};  
            case 'f0'
                f0 = varargin{i + 1};
            case 'cutEdges'
                cutEdges = varargin{i + 1};
            case 'd'
                d = varargin{i + 1};
            case 'p'
                p = varargin{i + 1};
            otherwise
                disp('problem with the agrument')
        end
    end
end

N = length(sig); sig = sig(:); 
sig = detrend(sig);
T = N / fs;
th = 0.001; % amplitude threshold for wavelet to be non-zero
om0 = 2 * pi * f0; %just for convenience denote circular central frequency

if strcmpi(wavelet, 'morlet')
    norm = @(om) (om0 / om)^(1-p);
    fwt = @(xi)(exp(-(1/2) * (om0-xi).^2) - exp(-(1/2) * (om0^2 + xi.^2)));  
    len = @(om) 2 * (om0 / om) * sqrt(-2 * log(sqrt(2 * pi) * th));
elseif strcmpi(wavelet, 'modMorlet')
    % modified Morlet, eq. 19 from J. Jamsek, A. Stefanovska, and P. V. E. McClintock, PrE 76, 046221, 2007
    omMin = 2 * pi / T; omMax = pi * fs;
    a = @(om) 2.^(1.8 * (om - omMin) / (omMax - omMin));
    c = 3.9487 * pi^(-0.25);
    % I try to find the length of Morlet wavelet...
    A = @(om) 1/sqrt(2 * pi) * a(om) * c;
    len = @(om) 2 * (om0 / om) * sqrt(-d * log(th / A(om)));
    norm = @(om) 1 / len(om) * (om0 ./ om);
    fwt = @(xi) a(xi) .* c * sqrt(d/2) .* (exp(-(d/4) * (om0-xi).^2) - exp(-om0^2 * d / 4 - xi.^2 / 2));
else
    error('Invalid wavelet name, Morlet or ModMorlet');
end

% frequencies
nq = ceil((N + 1) / 2); 
ff = [(0 : nq - 1), -fliplr(1 : N - nq)] * fs / N; 
ff = ff(:);

fx = fft(sig, N); 
fx(ff <= 0) = 0; 

wt = zeros(1, N) * NaN;
% wavelet transform
n1 = 0; n2 = n1; 
if cutEdges
    n1 = floor(fs * 0.5 * len(2 * pi * fr)); 
    n2 = n1; 
end
if n1 + n2 < N
    freqswf = ff * om0 / (2 * pi * fr); % for the wavelet
    fw = conj(fwt(2 * pi * freqswf)); 
    nanid = find(isnan(fw) | ~isfinite(fw));
    if ~isempty(nanid) % to avoid NaNs due to numerics, e.g. sin(0)/0
        fw(nanid) = conj(fwt(2 * pi * freqswf(nanid) + 10^(-14)));
        nanid = isnan(fw) | ~isfinite(fw); fw(nanid) = 0;
    end
    conv = fx .* fw(:); % convolution in the frequency domain
    out = norm(2 * pi * fr) * ifft(conv, N);
    wt(1 + n1 : N - n2) = out(1 + n1 : N - n2); % calculate WT at each time
end
end
