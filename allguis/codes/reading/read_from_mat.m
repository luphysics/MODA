function M = read_from_mat
%Reads from the file and creates a MATLAB variable to be used by the program

[filename,pathname,filterindex] = uigetfile('*.mat'); %only allows csv format files to be read

name = fullfile(pathname,filename)
M = load(name);
