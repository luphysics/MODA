% update for MODA 19.07.2018
% code compatible with Dmytro Iatsenko's wt.m
% author: Aleksandra Pidde a.pidde@gmail.com, a.pidde@lancaster.ac.uk
% mostly taken from Dmytro Iatsenko's wt.m code

function wt = wtAtf2(sig, fs, fr, opt)
% function calculating wavelet transform of signal sig at frequency fr
%
% INPUT:
% sig1: 		signal
% fs: 			sampling frequency
% fr:			frequency for each the wavelet is caluculated
% opt:			structure of optimal parameters returned by wt.m
%
% OUTPUT:
% wt:			wavelet transform at frequency fr

p = 1;
N = length(sig); sig = sig(:);
PadMode = opt.Padding;
L = N; fmin = opt.fmin; fmax = opt.fmax;

opt.fwt = @(xi)exp(-(q^2/2)*(log(xi).^2));

% ======== Signal preprocessing: detrending, filtering and padding =========
% [dflag] determines to do detrending and filtering before or after padding

dflag = 0;
if ~iscell(PadMode)
    if ~ischar(PadMode) && ~isempty(PadMode(PadMode ~= 0)), dflag = 1; end
    if strcmpi(PadMode, 'predictive') && fmin < 5 * fs / L, dflag = 1; end
else
    if ~ischar(PadMode{1}) && ~isempty(PadMode{1}(PadMode{1} ~= 0)), dflag = 1; end
    if ~ischar(PadMode{2}) && ~isempty(PadMode{2}(PadMode{2} ~= 0)), dflag = 1; end
    if strcmpi(PadMode{1}, 'predictive') && fmin < 5 * fs / L, dflag = 1; end
    if strcmpi(PadMode{2}, 'predictive') && fmin < 5 * fs / L, dflag = 1; end
end

% Detrend (subtract third-order polynomial fit) and filter first for usual padding

if strcmpi(opt.Preprocess, 'on') && dflag == 0
    [sig, fx, ff] = preprocess(sig, L, fs, fmin, fmax);
end
padleft = opt.PadLR{1}; padright = opt.PadLR{2};
sig = [padleft; sig; padright];
NL = length(sig);

%Detrend (subtract third-order polynomial fit) after padding for special cases
if strcmpi(opt.Preprocess, 'on') && dflag == 1
    [sig, fx, ff] = preprocess(sig, NL, fs, 0, fs / 2);
end

%Filtering of the padded signal
Nq = ceil((NL + 1) / 2);
ff = [(0 : Nq - 1), -fliplr(1 : NL - Nq)] * fs / NL; ff = ff(:); % frequencies in Fourier transform
fx = fft(sig, NL); fx(ff <= 0) = 0; % Fourier transform of a signal (set to zero at negative frequencies)
if strcmpi(opt.Preprocess,'on')
    fx(ff <= max([fmin, fs / L]) | ff >= fmax) = 0; % filter signal in a chosen frequency domain
end

coib1 = ceil(abs(opt.t1e * fs * (opt.ompeak./(2 * pi * fr))));
coib2 = ceil(abs(opt.t2e * fs * (opt.ompeak./(2 * pi * fr)))); %cone of influence edges

if (opt.t2e - opt.t1e) * opt.ompeak / (2 * pi * fmax) > L / fs
    coib1 = 0; coib2 = 0;
end

if coib1 == 0 && coib2 == 0
    n1 = floor((NL - L) / 2); n2 = ceil((NL - L) / 2);
else
    n1 = floor((NL - L) * coib1 / (coib1 + coib2));
    n2 = ceil((NL - L) * coib2 / (coib1 + coib2));
end
%plot(sig, '--');
% wavelet transform
wt = zeros(1, N) * NaN;

freqwf = ff * opt.ompeak / (2 * pi * fr); %frequencies for the wavelet
ii = find(freqwf > opt.xi1 / (2 * pi) & freqwf < opt.xi2 / (2 * pi)); %take into account only frequencies within the wavelet support
if ~isempty(opt.fwt)
    fw = conj(opt.fwt(2 * pi * freqwf(ii)));
    nid = find(isnan(fw) | ~isfinite(fw));
    if ~isempty(nid) %to avoid NaNs due to numerics, e.g. sin(0)/0
        fw(nid) = conj(opt.fwt(2 * pi * freqwf(ii(nid)) + 1e-14));
        nid = isnan(fw) | ~isfinite(fw); fw(nid) = 0;
    end
else
    twav = (1 / fs) * [-(1 : ceil((NL - 1) / 2)) + 1, NL + 1 -(ceil((NL - 1) / 2) + 1 : NL)]';
    timewf = (2 * pi * fr / opt.ompeak) * twav;
    jj = find(timewf > opt.t1 & timewf < opt.t2);
    tw = zeros(NL, 1); %take into account only times within the wavelet support
    tw(jj) = conj(opt.twf(timewf(jj)));
    nid = find(isnan(tw) | ~isfinite(tw));
    if ~isempty(nid) %to avoid NaNs due to numerics, e.g. sin(0)/0
        tw(nid) = conj(opt.twf(timewf(nid) + 1e-14));
        nid = isnan(tw) | ~isfinite(tw); tw(nid) = 0;
    end
    fw = (1 / fs) * fft(tw); fw = fw(ii);
end
cc = zeros(NL, 1); cc(ii) = fx(ii) .* fw(:); %convolution in the frequency domain
out = ((opt.ompeak / (2 * pi * fr))^(1-p)) * ifft(cc, NL); % calculate WT at each time
wt(1 : L) = out(1 + n1 : NL - n2);

if strcmpi(opt.CutEdges,'on')
    wt(1 : coib1) = nan;
    wt(end - coib2 : end) = nan;
end
end

function [newSig, fx, ff] = preprocess(sig, N, fs, fmin, fmax)
% Detrending
X = (1 : length(sig))' / fs; XM = ones(length(X), 4);
for pn = 1 : 3
    CX = X .^ pn;
    XM(:, pn + 1) = (CX - mean(CX)) / std(CX);
end
w = warning('off', 'all'); sig = sig - XM * (pinv(XM) * sig); warning(w);

% Filtering
fx = fft(sig, N); % Fourier transform of a signal
Nq = ceil((N + 1) / 2);
ff = [(0 : Nq - 1), -fliplr(1 : N - Nq)] * fs / N; ff = ff(:); % frequencies in Fourier transform
fx(abs(ff) <= max([fmin, fs / N]) | abs(ff) >= fmax) = 0; % filter signal in a chosen frequency domain
newSig = ifft(fx);
end
