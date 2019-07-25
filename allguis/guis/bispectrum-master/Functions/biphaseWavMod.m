
function [biamp, biphase] = biphaseWavMod(sig1, sig2, fs, f1, f2, varargin)
% Author:       Aleksandra Pidde 
%               a.pidde@lancaster.ac.uk, aleksandra.pidde@gmail.com

% date:         31.12.2017
% update:       8.01.2018
%
% calculating wavelet biphase and biamplitude
% sig1, sig2:   signals
% fs:           sampling frequency
% f1, f2:       pair of frequency to calculate biamplitude and biphase
%
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
%               params from eq 3.6), default 1
% OUTPUT
% biamp:        time-depandant biamplitude
% biphase:      time-depandant biphase

wavelet = 'morlet';
f0 = 1; 
cutEdges = true;
d = 2; 
p = 1;
if nargin >= 5 + 2
    for i = 1 : 2 : nargin - 5
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

auto = false;
if compareMatrix(sig1, sig2)
    auto = true;    
end
WT1 = wtAtfMod(sig1, fs, f1, 'f0', f0, 'wavelet', wavelet, 'cutEdges', cutEdges, 'p', p, 'd', d);
if auto && f1 == f2 
    WT2 = WT1;
else
    WT2 = wtAtfMod(sig2, fs, f2, 'f0', f0, 'wavelet', wavelet, 'cutEdges', cutEdges, 'p', p, 'd', d);
end
f3 = f1 + f2;
WT3 = wtAtfMod(sig2, fs, f3, 'f0', f0, 'wavelet', wavelet, 'cutEdges', cutEdges, 'p', p, 'd', d);

xx = WT1 .* WT2 .* conj(WT3);
biamp = abs(xx);
biphase = unwrap(angle(xx));
%biphase = angle(xx);
end
