
function [WT, freqs] = myWt(sig, fs, varargin)
% author: Aleksandra Pidde
%           a.pidde@lancaster.ac.uk, aleksandra.pidde@gmail.com
% date:     30.12.2017
% update:   8.01.2018
%
% calculating wavelets 
% INPUT
% sig:      signal
% fs:       sampling frequency [Hz]
%
% varargin (optional): 
% 'wavelet':type of wavelet, {'morlet', 'modMorlet'}, default 'morlet'
% 'f0':     tradeoff between the time and frequency resolutions, the higher, the
%           better frequency resolution, the lower the better time
%           resolution, default 1
% 'cutEdges':cutting of the cone of influence, {true, false}, default true
% 'fmin':   minimal frequency [Hz], bigger than 1 / T, where T is T of the
%           signal, default is such to fit one length of the wavelet in
%           small scale
% 'fmax':   maximal frequency [Hz], smaller than Nyquist freguency: fs / 2,
%           default is fs / 2
% 'nv':     number of voices, frequency discretization, meaning the next 
%           frequency equals previous one multiplied on [2^(1/nv)], default
%           16 
% 'd':      Gaussian decay in Morlet wavelet, according to
%           correction, default 2
% 'p':      normalization factor (Kaiser, A friendly guide to wavelets,
%           params from eq 3.6), default 1/2, preserving the energy
%
% OUTPUT
% WT:       wavelet transform
% freq:     frequency vector

N = length(sig); sig = sig(:); 
sig = detrend(sig);
T = N / fs;
th = 0.001; % amplitude threshold for wavelet to be non-zero
% default
wavelet = 'morlet';
f0 = 1; 
cutEdges = true;
fmin = {};
fmax = fs / 2; 
nv = 16;
p = 0.5; 
d = 2;
% if defined by user
if nargin >= 2 + 2
    for i = 1 : 2 : nargin - 2
        switch varargin{i}
            case 'wavelet'
                wavelet = varargin{i + 1};  
            case 'f0'
                f0 = varargin{i + 1};
            case 'cutEdges'
                cutEdges = varargin{i + 1};
            case 'fmin'
                fmin = varargin{i + 1};
            case 'fmax'
                fmax = varargin{i + 1};
            case 'nv'
                nv = varargin{i + 1};
            case 'd'
                d = varargin{i + 1};
            case 'p'
                p = varargin{i + 1};
            otherwise
                disp('problem with the agrument')
        end
    end
end

om0 = 2 * pi * f0; 
if strcmpi(wavelet, 'morlet')
    norm = @(om) (om0 / om)^(1-p);
    fwt = @(xi)(exp(-(1/2) * (om0-xi).^2) - exp(-(1/2) * (om0^2 + xi.^2)));  
    len = @(om) 2 * (om0 / om) * sqrt(-2 * log(sqrt(2 * pi) * th));
    if isempty(fmin)
        fmin = 2 * f0 / T * sqrt(2) * sqrt(-log(sqrt(2 * pi) * th));
    end
elseif strcmpi(wavelet, 'modMorlet')
    % modified Morlet, eq. 19 from J. Jamsek, A. Stefanovska, and P. V. E. McClintock, PrE 76, 046221, 2007
    omMin = 2 * pi / T; omMax = 2 * pi * fmax;
    a = @(om) 2.^(1.8 * (om - omMin) / (omMax - omMin));
    c = 3.9487 * pi^(-0.25);
    % I try to find the length of Morlet wavelet...
    A = @(om) 1/sqrt(2 * pi) * a(om) * c;
    len = @(om) 2 * (om0 / om) * sqrt(-d * log(th / A(om)));
    norm = @(om) 1 / len(om) * (om0 ./ om);
    fwt = @(xi) a(xi) .* c * sqrt(d/2) .* (exp(-(d/4) * (om0-xi).^2) - exp(-om0^2 * d / 4 - xi.^2 / 2));
    if isempty(fmin)
        fmin = 2 * f0 / T * sqrt(2) * sqrt(-log(sqrt(2 * pi) * th));
    end
else
    error('Invalid wavelet name, Morlet or ModMorlet');
end

% frequencies
freqs = 2.^((ceil(nv * log2(fmin)) : floor(nv * log2(fmax)))' / nv);
nf = length(freqs); % number of frequencies
%freqs = linspace(fmin, fmax, nf); %checkpoint what happens when generate for linear axis
nq = ceil((N + 1) / 2); 
ff = [(0 : nq - 1), -fliplr(1 : N - nq)] * fs / N; 
ff = ff(:);

fx = fft(sig, N); 
fx(ff <= 0) = 0; 

%%
nor = zeros(1, nf);
% wavelet transform
WT = zeros(nf, N) * NaN; 
for sn = 1 : nf
    %cone of influence
    n1 = 0; n2 = n1; 
    if cutEdges
        n1 = floor(fs * 0.5 * len(2 * pi * freqs(sn))); 
        n2 = n1; 
    end
    if n1 + n2 < N
        freqswf = ff * om0 / (2 * pi * freqs(sn)); % for the wavelet
        fw = conj(fwt(2 * pi * freqswf)); 
        nanid = find(isnan(fw) | ~isfinite(fw));
        if ~isempty(nanid) % to avoid NaNs due to numerics, e.g. sin(0)/0
            fw(nanid) = conj(fwt(2 * pi * freqswf(nanid) + 10^(-14)));
            nanid = isnan(fw) | ~isfinite(fw); fw(nanid) = 0;
        end
        conv = fx .* fw(:); % convolution in the frequency domain
        out = norm(2 * pi * freqs(sn)) * ifft(conv, N); % calculate WT at each time
        %cone of influence
        WT(sn, 1 + n1 : N - n2) = out(1 + n1 : N - n2);
    end
end

end
