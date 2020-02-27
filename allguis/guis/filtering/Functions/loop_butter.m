%function looping butterworth filter till reaching the maximum order
% author: Valentina Ticcinelli
%output:
%s_out= filtered signal
%or= optimal order

%input:
%s_in= signal to be filtered
%band=[f_min f_max] passing band
%fs= sample frequency

function [s_out,or]=loop_butter(s_in,band,fs)

or=1;
s_out=bandpass_butter(s_in,or,band(1),band(2), fs);

max_out=max(abs(s_out));

while max(abs(s_out))<10*max_out
    or=or+1;
s_out=bandpass_butter(s_in,or,band(1),band(2), fs);
end
or=or-2;

if or<1
    or=1;
else
end
s_out=bandpass_butter(s_in,or,band(1),band(2), fs);


function [d]=bandpass_butter(c,n,flp,fhi,fs) 

% [d]=bandpass(c,flp) 
% 
% bandpass a time series with a nth order butterworth filter 
% 
% c = input time series 
% n= order for the filter
% flp = lowpass corner frequency of filter 
% fhi = hipass corner frequency 
% fs = sampling frequency
 

fnq=1/2*fs;  % Nyquist frequency 
Wn=[flp/fnq fhi/fnq];    % butterworth bandpass non-dimensional frequency 
[b,a]=butter(n,Wn); % construct the filter 
d=filtfilt(b,a,c); % zero phase filter the data 


