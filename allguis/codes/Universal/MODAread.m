% MODA data loading function

function [handles,sig,E]=MODAread(handles,type,varargin)

% Parse varargin for whether the number of signals should be even.
% This is used in phase coherence, bispectrum analysis and 
% Bayesian inference.
even = false;
n = length(varargin);
for k = 1:n
    if strcmp("even", varargin{k})
        even = true;
        break;
    end
end

E=1;
set(handles.status,'String','Importing Signal...');  % Update status

[filename,pathname] = uigetfile('*.*');  % Request file location from user
name = fullfile(pathname,filename);

if filename==0
    sig=0;
    return;
end

try
    sig = load(name); % Load data.
catch exception
    % Catch exception when opening Excel CSV files with BOM.
    sig = readmatrix(name);
end

if isstruct(sig) % If loaded signal is a MATLAB structure, convert to array
    sig=struct2array(sig);
else
end

handles.sampling_freq = str2double(cell2mat(newid(['Enter the sampling frequency of the data (',filename,') in Hz'])));
fs = handles.sampling_freq;

if isnan(fs)
    errordlg('Sampling frequency must be specified')
    E=0;
    return;
end

choice = questdlg('Select Orientation of Data set?', ...
    'Data Import','Column wise','Row wise','default');
switch choice
    case 'Column wise'
        sig = sig';
        
end

if isempty(choice)
    errordlg('Data set orientation must be specified')
    E=0;
    return;
end

num_signals = length(sig(:,1));

% If there are an odd number of signals but an even number must be 
% supplied, remove the last one.
if even && mod(num_signals, 2) ~= 0 
    sig = sig(1:end-1,:);
end

handles.sig = sig;
handles.sig_cut=handles.sig;
handles.sig_pp=sig;
time = linspace(0,length(handles.sig)/handles.sampling_freq,length(handles.sig));
handles.time_axis = time;
handles.time_axis_cut = time;
handles.xl=[handles.time_axis_cut(1) handles.time_axis_cut(end)];

if type==1
    N=size(sig);
    if N(1)==1;
        handles.sig=[handles.sig;handles.sig];
        handles.sig_cut=[handles.sig;handles.sig];
        handles.sig_pp=[handles.sig;handles.sig];
    else
        
        %% Plot time series
        linkaxes([handles.time_series_1 handles.time_series_2],'x'); % Ensures axis limits are identical for both plots
        plot(handles.time_series_1,handles.time_axis,handles.sig(1,:),'color',handles.linecol(1,:));
        xlim(handles.time_series_1,[0 handles.time_axis(end)]);
        plot(handles.time_series_2,handles.time_axis,handles.sig(1+size(handles.sig,1)/2,:),'color',handles.linecol(1,:));
        xlim(handles.time_series_2,[0 handles.time_axis(end)]);
        xlabel(handles.time_series_2,'Time (s)');
        ylabel(handles.time_series_1,'Sig 1');
        ylabel(handles.time_series_2,'Sig 2');
        
        if   mod(N(1),2)==1;
            errordlg('Number of data sets must be even','Data Error');
            E=0;
            return;
            
        end
        
        %% Create signal list
        if isfield(handles,'signal_list')
            list = cell(size(sig,1)/2,1);
            list{1,1} = 'Signal Pair 1';
            
            for i = 2:size(sig,1)/2
                list{i,1} = sprintf('Signal Pair %d',i);
            end
            
            set(handles.signal_list,'String',list);
        else
        end
    end
else
    
    %% Plot time series
    plot(handles.time_series,handles.time_axis,sig(1,:),'color',handles.linecol(1,:));
    xlim(handles.time_series,[0 handles.time_axis(end)]);
    xlabel(handles.time_series,'Time (s)');
    ylabel(handles.time_series,'Sig');
    
end




set(handles.plot_TS,'Enable','on')

set(handles.status,'String','Select data and define parameters');


