% Copy of biphaseWavNew to be used as a Python library.

function [biamp, biphase] = biphaseWavPython(sig1, sig2, fs, f0, fr, opt)
% calculating wavelet biamplitude and biphase of two signal
%
% INPUT:
% sig1, sig2:   signals
% fs:           sampling frequency
% fr:           pair of frequencies, [f1, f2], f1 is from sig1, f2 from sig2
% opt:          structure of optimal parameters returned by wt.m
%
% OUTPUT:
% biamp:        wavelet biamplitude
% biphase:      wavelet biphase

q = 2*pi*f0;
opt.fwt = @(xi)exp(-(q^2/2)*(log(xi).^2));

opt.PadLR = {};
opt.PadLR{1} = opt.PadLR1;
opt.PadLR{2} = opt.PadLR2;

opt.twf = {};
opt.twf{1} = complex(opt.twf1r, opt.twf1i);
opt.twf{2} = complex(opt.twf2r, opt.twf2i);

f1 = fr(1);
f2 = fr(2);
f3 = f1 + f2;

wt1 = wtAtf2Python(sig1, fs, f1, opt);
wt2 = wtAtf2Python(sig2, fs, f2, opt);
wt3 = wtAtf2Python(sig2, fs, f3, opt);

xx = wt1 .* wt2 .* conj(wt3);
biamp = abs(xx);
biphase = unwrap(angle(xx));
end