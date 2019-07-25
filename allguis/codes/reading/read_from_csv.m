function M = read_from_csv
%Reads from the file and creates a MATLAB variable to be used by the program

[filename,pathname,filterindex] = uigetfile('*.csv'); %only allows csv format files to be read
name = fullfile(pathname,filename)

M = csvread(name);
