
function [Bisp, freqs, Norm, WT1, WT2] = bispecWavMod(sig1, sig2, fs,varargin)
% Author:       Aleksandra Pidde 
%               a.pidde@lancaster.ac.uk, aleksandra.pidde@gmail.com
% date:         30.12.2017
% update:       8.01.2018
%
% calculating wavelet bispectrum
% sig1, sig2:   signals
% fs:           sampling frequency
%
% varargin (optional): 
% 'wavelet':    type of wavelet, {'morlet', 'modMorlet'}, default 'morlet'
% 'f0':         tradeoff between the time and frequency resolutions, the higher, the
%               better frequency resolution, the lower the better time resolution,
%               default 1 
% 'fmin':       minimal frequency [Hz], bigger than 1 / T, where T is T of
%               the signal, default is such to fit one length of the wavelet in
%               big scale
% 'fmax':       maximal frequency [Hz], smaller than Nyquist freguency: 
%               fs / 2, default is fs / 2
% 'nv':         number of voices, frequency discretization, meaning the next 
%               frequency equals previous one multiplied on [2^(1/nv)],
%               default 16
% 'cutEdges':   cutting of the cone of influence, {true, false}, default
%               true
% 'd':          Gaussian decay in Morlet wavelet, according to
%               correction, default 2
% 'p':          normalization factor (Kaiser, A friendly guide to wavelets,
%               params from eq 3.6), default 1
%
% OUTPUT
% Bisp:         bispectrum
% freq:         vector of frequencies
% Norm:         normalisation for bicoherence
try

bstype={'111';'222';'122';'211'};
N = length(sig1); sig1 = sig1(:); 
sig1 = detrend(sig1);
sig2 = sig2(:); sig2 = detrend(sig2);
T = N / fs;
th = 0.001; % amplitude threshold for wavelet to be non-zero
% default
wavelet = 'morlet';
f0 = 1; 
cutEdges = true;
fmin = {};
fmax = fs / 2; 
nv = 16;
p = 1; 
d = 2;
wbar=0;
% if defined by user
if nargin >= 2 + 3
    for i = 1 : 2 : nargin - 3
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
            case 'handles'
                handles=varargin{i+1};
            case 'hObject'
                hObject=varargin{i+1};
            case 'num'
                numb=varargin{i+1};
            case 'wbar'
                wbar=varargin{i+1};
            otherwise
                disp('problem with the argument')
        end
    end
end

if wbar==1
handles.h = waitbar(0,'Calculating bispectrum...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(handles.h,'canceling',0)
    guidata(hObject,handles); 
else
end

    
om0 = 2 * pi * f0; %just for convenience denote circular central frequency
if isempty(fmin)
    fmin = 2 * f0 / T * sqrt(2) * sqrt(-log(sqrt(2 * pi) * th));
end

if isnan(fmin)
    fmin = 2 * f0 / T * sqrt(2) * sqrt(-log(sqrt(2 * pi) * th));
end
if strcmpi(wavelet, 'morlet')
    norm = @(om) (om0 / om)^(1-p);
    fwt = @(xi)(exp(-(1/2) * (om0-xi).^2) - exp(-(1/2) * (om0^2 + xi.^2)));  
    len = @(om) 2 * (om0 / om) * sqrt(-2 * log(sqrt(2 * pi) * th));
    
elseif strcmpi(wavelet, 'modMorlet')
    % modified Morlet, eq. 19 from J. Jamsek, A. Stefanovska, and P. V. E. McClintock, PrE 76, 046221, 2007
    omMin = 2 * pi * fmin; omMax = 2 * pi * fmax;
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

auto = false;
if compareMatrix(sig1, sig2)
    auto = true;    
end
[WT1, freqs] = myWt(sig1, fs, 'f0', f0, 'fmin', fmin, 'fmax', fmax, 'wavelet', wavelet, 'nv', nv, 'cutEdges', cutEdges, 'p', p, 'd', d);
if auto 
    WT2 = WT1;
else
    [WT2, freqs] = myWt(sig2, fs, 'f0', f0, 'fmin', fmin, 'fmax', fmax, 'wavelet', wavelet, 'nv', nv, 'cutEdges', cutEdges, 'p', p, 'd', d);
end

nf = length(freqs); %number of frequencies
Bisp = NaN * zeros(nf, nf);
Norm = ones(nf, nf);

nq = ceil((N + 1) / 2); 
ff = [(0 : nq - 1), -fliplr(1 : N - nq)] * fs / N; 
ff = ff(:);

fx = fft(sig2, N); 
fx(ff <= 0) = 0; 

for j = 1 : nf
    if wbar==1
    if getappdata(handles.h,'canceling')
        guidata(hObject,handles);
                break;
    end
    else
    end

    kstart = 1;
    if auto
        kstart = j;
    end
    for k = kstart : nf
        if wbar==1
        if getappdata(handles.h,'canceling')
            guidata(hObject,handles);
                break;
        end
        else
        end

        f3 = freqs(j) + freqs(k);
        bigger = max([j k]);
        idx3 = find(freqs >= f3, 1);
        if (f3 <= freqs(end)) && (freqs(idx3 - 1) > freqs(bigger)) 
            WTdat = wtAtfMod(sig2, fs, f3, 'f0', f0, 'wavelet', wavelet, 'cutEdges', cutEdges, 'p', p, 'd', d);
            WTdat = WTdat(:).'; % make sure it is horizontal vector
            xx = WT1(j, :) .* WT2(k, :) .* conj(WTdat);
            ss = nanmean(xx);
            Bisp(j, k) = ss;
            
            norm1 = nanmean(abs(WT1(j, :) .* WT2(k, :)).^2);
            norm2 = nanmean(abs(WTdat).^2);
            a = sqrt(norm1 * norm2);
            Norm(j, k) = a;
        end
    end
    if wbar==1
    waitbar(j / length(freqs),handles.h,sprintf(['Calculating Bispectrum ',bstype{numb},' (%d/%d)'],j,length(freqs)));
    else
    end
end
if wbar==1
delete(handles.h); 
else
end



catch e
    errordlg(e.message,'Error'); 
    if wbar==1
    delete(handles.h); 
    else
    end
    rethrow(e)
end
end
