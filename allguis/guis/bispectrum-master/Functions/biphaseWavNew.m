% update for MODA: 19.07.2018
% code compatible with Dmytro Iatsenko's wt.m
% author: Aleksandra Pidde aleksandra.pidde@gmail.com,
% a.pidde@lancaster.ac.uk

function [biamp, biphase] = biphaseWavNew(sig1, sig2, fs, fr, opt)
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

f1 = fr(1);
f2 = fr(2);
f3 = f1 + f2;
wt1 = wtAtf2(sig1, fs, f1, opt);
wt2 = wtAtf2(sig2, fs, f2, opt);
wt3 = wtAtf2(sig2, fs, f3, opt);

xx = wt1 .* wt2 .* conj(wt3);
biamp = abs(xx);
biphase = unwrap(angle(xx));
end


