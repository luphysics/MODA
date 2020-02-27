function [d] = bandpass_butter(c,or,flp,fhi,fs) 
% 
% [d]=bandpass(c,flp) 
% 
% bandpass a time series with a 2nd order butterworth filter 
% 
% c = input time series 
% or-order for the filter
% flp = lowpass corner frequency of filter 
% fhi = hipass corner frequency 
% fs = sampling frequency
% 
n=or;      % 2nd order butterworth filter 
fnq=1/2*fs;  % Nyquist frequency 
Wn=[flp/fnq fhi/fnq];    % butterworth bandpass non-dimensional frequency 
[b,a]=butter(n,Wn); % construct the filter 
d=filtfilt(b,a,c); % zero phase filter the data 
return;