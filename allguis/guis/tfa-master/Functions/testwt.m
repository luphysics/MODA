fs = 20;
t = 0:1/fs:50-1/fs;

signal = cos(2*pi*3*t + .75 * sin(2*pi*t/5));

result = wt(signal, fs);