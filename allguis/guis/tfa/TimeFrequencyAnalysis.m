%Version 1.01
%********************************************************************************
%*************************** Time-Frequency Analysis GUI ************************
%********************************************************************************
%---------------------------Credits---------------------------------------
% Wavelet and windowed Fourier transform: Dmytro Iatsenko
%
%----------------------------Documentation--------------------------------
% Reads a single or matrix of signals in any format readable by MATLAB.
% User can select the part of the signal they want to use, and calculate the wavelet
% transform or windowed Fourier transform of that part.
% Displays the amplitude/power surface plot and the time-averaged plot.
% Also contains save options for the graphs and data from the wavelet and
% windowed Fourier transforms.



function varargout = TimeFrequencyAnalysis(varargin)
%      TIMEFREQUENCYANALYSIS MATLAB code for TimeFrequencyAnalysis.fig
%      TIMEFREQUENCYANALYSIS, by itself, creates a new TIMEFREQUENCYANALYSIS or raises the existing
%      singleton*.
%
%      H = TIMEFREQUENCYANALYSIS returns the handle to a new TIMEFREQUENCYANALYSIS or the handle to
%      the existing singleton*.
%
%      TIMEFREQUENCYANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIMEFREQUENCYANALYSIS.M with the given input arguments.
%
%      TIMEFREQUENCYANALYSIS('Property','Value',...) creates a new TIMEFREQUENCYANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TimeFrequencyAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TimeFrequencyAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help TimeFrequencyAnalysis

% Last Modified by GUIDE v2.5 26-Mar-2018 14:35:48
%*************************************************************************%
%                BEGIN initialization code - DO NOT EDIT                  %
%                ----------------------------------------                 %
%*************************************************************************%
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TimeFrequencyAnalysis_OpeningFcn, ...
    'gui_OutputFcn',  @TimeFrequencyAnalysis_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%*************************************************************************%
%                END initialization code - DO NOT EDIT                    %
%*************************************************************************%


function TimeFrequencyAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles=MODAsettings(hObject, handles);

% Disable plot functions on startup
set(handles.plot_TS,'Enable','off')
set(handles.save_3dplot,'Enable','off')
set(handles.save_both_plot,'Enable','off')
set(handles.save_avg_plot,'Enable','off')
set(handles.save_mm_plot,'Enable','off')
set(handles.kaisera,'Enable','off')

% Disable save functions on startup
set(handles.mat_save,'Enable','off')
set(handles.csv_save,'Enable','off')
set(handles.save_WT_coeff,'Enable','off')
set(handles.save_session,'Enable','off')

drawnow;
handles.output = hObject;
guidata(hObject, handles);

function file_read_Callback(hObject, eventdata, handles)
% Loads user data when selected from File --> Load time series

[handles,A]=MODAreadcheck(handles);
if A==1
    
    
    %------- Resetting axes and clearing data from previous ---
    %------- calculations. ------------------------------------
    
    
    cla(handles.cum_avg,'reset');
    cla(handles.plot3d,'reset');
    cla(handles.plot_pow,'reset');
    set(handles.cum_avg,'visible','off');
    set(handles.plot3d,'visible','on');
    set(handles.plot_pow,'visible','on');
    uistack(handles.plot3d,'top');
    uistack(handles.plot_pow,'top');
    cla(handles.time_series,'reset')
    cla(handles.plot_pp,'reset')
    set(handles.signal_list,'String',[]);
    
    if isfield(handles, 'freqarr');handles = rmfield(handles, 'freqarr');else end
    if isfield(handles, 'sig');handles = rmfield(handles, 'sig');else end
    if isfield(handles, 'sig_cut');handles = rmfield(handles, 'sig_cut');else end
    if isfield(handles, 'sig_pp');handles = rmfield(handles, 'sig_pp');else end
    if isfield(handles, 'time_axis'); handles = rmfield(handles, 'time_axis');else end
    if isfield(handles, 'time_axis_us');handles = rmfield(handles, 'time_axis_us');else end
    if isfield(handles, 'pow_arr');handles = rmfield(handles, 'pow_arr');else end
    if isfield(handles, 'amp_arr');handles = rmfield(handles, 'amp_arr');else end
    if isfield(handles, 'pow_WT');handles = rmfield(handles, 'pow_WT');else end
    if isfield(handles, 'amp_WT');handles = rmfield(handles, 'amp_WT');else end
    if isfield(handles, 'WT');handles = rmfield(handles, 'WT');else end
    if isfield(handles, 'wopt');handles = rmfield(handles, 'wopt');else end
    if isfield(handles, 'sampling_freq');handles = rmfield(handles, 'sampling_freq');else end
    if isfield(handles, 'peak_value');handles = rmfield(handles, 'peak_value');else end
    
    % Load data
    [handles,sig,E]=MODAread(handles,0);
    
    if sig==0
    elseif E==0
    else
        % Create signal list
        list = cell(size(sig,1)+1,1);
        list{1,1} = 'Signal 1';
        for i = 2:size(sig,1)
            list{i,1} = sprintf('Signal %d',i);
        end
        set(handles.signal_list,'String',list);
        list{size(sig,1)+1,1} = sprintf('Average Plot (All)');
        set(handles.signal_list,'String',list);
        
        globalfontsize = 12; % Do not edit this line manually. See scripts/fontsize.py.
        
        % Plot time series
        plot(handles.time_series,handles.time_axis,sig(1,:),'color',handles.linecol(1,:));
        set(handles.time_series, "FontSize", globalfontsize);

        xlim(handles.time_series,[0,size(sig,2)./handles.sampling_freq]);
        xlabel(handles.time_series,'Time (s)');
        guidata(hObject,handles);
        refresh_limits_Callback(hObject, eventdata, handles);
        cla(handles.plot_pp,'reset');
        preprocess_Callback(hObject, eventdata, handles);
        xlabel(handles.time_series,'Time (s)');
        set(handles.status,'String','Select data and continue with transform');
        set(handles.signal_length,'String',strcat(num2str(size(sig,2)/handles.sampling_freq/60),' minutes'));
    end
else
    return;
end


function varargout = TimeFrequencyAnalysis_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function preprocess_Callback(hObject, eventdata, handles)

% Get the x-limits.
lim = csv_to_mvar(get(handles.xlim,"String"));
times = handles.time_axis;

% Indices between x-limits.
indices =  times >= lim(1) & times <= lim(2);

cla(handles.plot_pp,'reset');
L = size(handles.sig,2);
signal_selected = get(handles.signal_list, 'Value');
fs = handles.sampling_freq;
fmax = str2double(get(handles.max_freq,'String'));
fmin = str2double(get(handles.min_freq,'String'));

% Modify signal to be in correct range; make backup to restore later.
sigBackup = zeros(size(handles.sig));
sigBackup(:) = handles.sig;
handles.sig = handles.sig(:,indices);

handles.sig_pp = cell(size(handles.sig,1),1);

for i = 1:size(handles.sig,1)
    sig = handles.sig(i,:);
    
    % Detrending
    X=(1:length(sig))'/fs; XM=ones(length(X),4);
    
    for pn=1:3
        CX=X.^pn;
        XM(:,pn+1)=(CX-mean(CX))/std(CX);
    end
    sig = sig(:);
    w=warning('off','all');
    new_signal=sig-XM*(pinv(XM)*sig);
    warning(w);
    
    % Filtering
    fx=fft(new_signal,L);
    
    Nq=ceil((L+1)/2);
    ff=[(0:Nq-1),-fliplr(1:L-Nq)]*fs/L;
    ff=ff(:);
    
    fx(abs(ff)<=max([fmin,fs/L]) | abs(ff)>=fmax)=0;
    
    result = ifft(fx)';
    handles.sig_pp{i,1} = result;
    
end

handles.sig = sigBackup; % Restore full signal.

%%%% Plotting %%%%
globalfontsize = 12; % Do not edit this line manually. See scripts/fontsize.py.

% Plot original signal on preprocessing plot.
plot(handles.plot_pp,handles.time_axis,handles.sig(signal_selected,:),'color',handles.linecol(1,:));
hold(handles.plot_pp,'on');

% Plot preprocessed signal on preprocessing plot.
pp_times = times(indices);
pp_sig = handles.sig_pp{signal_selected,1};
pp_sig = pp_sig(1:size(find(1 == indices),2));
plot(handles.plot_pp, pp_times, pp_sig,'color',handles.linecol(2,:));

legend(handles.plot_pp,{'Original','Pre-Processed'},'FontSize',globalfontsize-2,'Location','Best','units','normalized');
xlim(handles.plot_pp,[0,size(handles.sig,2)./fs]);
set(handles.plot_pp,'FontUnits','points','FontSize',globalfontsize);
xlabel(handles.plot_pp,{'Time (s)'},'FontUnits','points','FontSize',globalfontsize-2);
xl = csv_to_mvar(get(handles.xlim,'String'));
xl = xl.*fs;
xl(2) = min(xl(2),size(handles.sig,2));
xl(1) = max(xl(1),1);
xl = xl./fs;
set(handles.plot_pp,'xlim',[xl(1) xl(2)]);
guidata(hObject,handles);
drawnow;

function signal_list_Callback(hObject, eventdata, handles)
% Selecting signal and calling other necessary functions
%handles.p=[];
signal_selected = get(handles.signal_list, 'Value');

if any(signal_selected == size(handles.sig,1)+1)
    set(handles.signal_list,'Max',size(handles.sig,1));
    set(handles.save_WT_coeff,'Enable','off')
else
    if size(signal_selected,2) == 1
        set(handles.signal_list,'Max',1);
    else
        set(handles.signal_list, 'Value', 1);
        set(handles.signal_list,'Max',1);
        drawnow;
        xyplot_Callback(hObject, eventdata, handles);
    end
end

if any(signal_selected ~= size(handles.sig,1)+1) && length(signal_selected) == 1
    set(handles.save_WT_coeff,'Enable','on')
    
    globalfontsize = 12; % Do not edit this line manually. See scripts/fontsize.py.

    plot(handles.time_series, handles.time_axis, handles.sig(signal_selected,:),'color',handles.linecol(1,:));
    set(handles.time_series, "FontSize", globalfontsize);

    xl = csv_to_mvar(get(handles.xlim, 'String'));
    xlim(handles.time_series, xl);
    xlabel(handles.time_series, 'Time (s)');
    refresh_limits_Callback(hObject, eventdata, handles);
    cla(handles.plot_pp, 'reset');
    preprocess_Callback(hObject, eventdata, handles);
    xlabel(handles.time_series, 'Time (s)');
    set(handles.status, 'String', 'Select data and continue with transform');
    
    if isfield(handles,'amp_WT')
        xyplot_Callback(hObject, eventdata, handles);
    end
    intervals_Callback(hObject, eventdata, handles)
elseif any(signal_selected == size(handles.sig,1)+1)
    xyplot_Callback(hObject, eventdata, handles);
    intervals_Callback(hObject, eventdata, handles)
end


function wavlet_transform_Callback(hObject, eventdata, handles)
handles.currsig=[];
handles=MODATFAcalc(hObject,eventdata,handles,1);
if ~handles.failed
    intervals_Callback(hObject, eventdata, handles)
    xyplot_Callback(hObject, eventdata, handles);
end
guidata(hObject,handles);

function wt_single_Callback(hObject, eventdata, handles)
handles.currsig=get(handles.signal_list,'Value');
handles=MODATFAcalc(hObject,eventdata,handles,2);
if ~handles.failed
    intervals_Callback(hObject, eventdata, handles)
    xyplot_Callback(hObject, eventdata, handles);
end


function xyplot_Callback(hObject, eventdata, handles)
% Plots all figures
signal_selected = get(handles.signal_list,'Value');
if any(signal_selected == size(handles.sig,1)+1) && isfield(handles,'freqarr')
    cla(handles.plot3d,'reset');
    cla(handles.plot_pow,'reset');
    cla(handles.cum_avg,'reset');
    cla(handles.time_series,'reset');
    set(handles.plot3d,'visible','off');
    set(handles.plot_pow,'visible','off');
    set(handles.cum_avg,'visible','on');
    set(handles.save_3dplot,'Enable','off')
    set(handles.save_both_plot,'Enable','off')
    set(handles.save_avg_plot,'Enable','off')
    set(handles.save_mm_plot,'Enable','on')

    globalfontsize = 12; % Do not edit this line manually. See scripts/fontsize.py.
    
    hold(handles.cum_avg,'on');
    uistack(handles.cum_avg,'top');
    size(handles.sig,1);
    if(handles.plot_type == 1)
        if size(handles.sig,1) == 1
            plot(handles.cum_avg, handles.freqarr,cell2mat(handles.pow_arr),'-','Linewidth',3,'color',handles.linecol(1,:)); % Linecolors edited 04/09/2017 - GL
            plot(handles.cum_avg, handles.freqarr,cell2mat(handles.pow_arr),'--','Linewidth',3,'color',handles.linecol(2,:));
        else
            plot(handles.cum_avg, handles.freqarr, mean(cell2mat(handles.pow_arr)),'-','Linewidth',3,'color',handles.linecol(1,:));
            plot(handles.cum_avg, handles.freqarr, median(cell2mat(handles.pow_arr)),'--','Linewidth',3,'color',handles.linecol(2,:));
        end
        ylabel(handles.cum_avg,'Average Power','FontUnits','points','FontSize',globalfontsize);
        xlabel(handles.cum_avg,'Frequency (Hz)','FontUnits','points','FontSize',globalfontsize);
    else
        if size(handles.sig,1) == 1
            plot(handles.cum_avg, handles.freqarr, cell2mat(handles.amp_arr),'-','Linewidth',3,'color',handles.linecol(1,:));
            plot(handles.cum_avg, handles.freqarr, cell2mat(handles.amp_arr),'--','Linewidth',3,'color',handles.linecol(2,:));
        else
            plot(handles.cum_avg, handles.freqarr, mean(cell2mat(handles.amp_arr)),'-','Linewidth',3,'color',handles.linecol(1,:));
            plot(handles.cum_avg, handles.freqarr, median(cell2mat(handles.amp_arr)),'--','Linewidth',3,'color',handles.linecol(2,:));
        end
        ylabel(handles.cum_avg,'Average Amplitude','FontUnits','points','FontSize',globalfontsize);
        xlabel(handles.cum_avg,'Frequency (Hz)','FontUnits','points','FontSize',globalfontsize);
    end
    handles.leg1={'Mean','Median'};
    % legend(handles.cum_avg,handles.leg1,'FontUnits','points','FontSize',10)
    legend(handles.cum_avg,handles.leg1,'FontSize',globalfontsize)
    ind=2;
    ls=1;
    sty='-';
    for i = 1:size(signal_selected,2)
        ind=ind+1;
        if ind>7
            ind=1;
            ls=ls*-1;
            if ls<0
                sty='-.';
            else
                sty='-';
            end
        else
        end
        if(handles.plot_type == 1 && signal_selected(i) <= size(handles.sig,1))
            plot(handles.cum_avg, handles.freqarr, handles.pow_arr{signal_selected(i),1},sty,'color',handles.linecol(ind,:),'LineWidth',handles.line2width);
            ylabel(handles.cum_avg,'Average Power','FontUnits','points','FontSize',globalfontsize);
            xlabel(handles.cum_avg,'Frequency (Hz)','FontUnits','points','FontSize',globalfontsize);
            [M,I] = max(handles.pow_arr{signal_selected(i),1});
            handles.leg1{i+2}=['Signal ',num2str(signal_selected(i))];
            legend(handles.cum_avg,handles.leg1)
            
        elseif signal_selected(i) <= size(handles.sig,1)
            plot(handles.cum_avg, handles.freqarr, handles.amp_arr{signal_selected(i),1},sty,'color',handles.linecol(ind,:),'LineWidth',handles.line2width);
            ylabel(handles.cum_avg,'Average Amplitude','FontUnits','points','FontSize',globalfontsize);
            xlabel(handles.cum_avg,'Frequency (Hz)','FontUnits','points','FontSize',globalfontsize);
            [M,I] = max(handles.amp_arr{signal_selected(i),1});
            handles.leg1{i+2}=['Signal ',num2str(signal_selected(i))];
            legend(handles.cum_avg,handles.leg1)
            
        end
        
    end
    grid(handles.cum_avg,'off');
    box(handles.cum_avg,'on');
    title(handles.cum_avg,'Transform average for all signals')
    if handles.calc_type == 1
        set(handles.cum_avg,'xscale','log');
    else
        set(handles.cum_avg,'xscale','linear');
    end
    
    xlim(handles.cum_avg,[min(handles.freqarr) max(handles.freqarr)]);
elseif isfield(handles,'freqarr')
    set(handles.save_3dplot,'Enable','on')
    set(handles.save_both_plot,'Enable','on')
    set(handles.save_avg_plot,'Enable','on')
    set(handles.save_mm_plot,'Enable','off')
    
    cla(handles.cum_avg,'reset');
    cla(handles.plot3d,'reset');
    cla(handles.plot_pow,'reset');
    set(handles.cum_avg,'visible','off');
    set(handles.plot3d,'visible','on');
    set(handles.plot_pow,'visible','on');
    uistack(handles.plot3d,'top');
    uistack(handles.plot_pow,'top');
    
    if(handles.plot_type == 1)
        WTpow = handles.pow_WT{signal_selected,1};
        handles.peak_value = max(WTpow(:))+.1;
        pcolor(handles.plot3d, handles.time_axis_us , handles.freqarr, WTpow);
        
        colorbar(handles.plot3d)
        plot(handles.plot_pow, handles.pow_arr{signal_selected,1}, handles.freqarr,'-k','LineWidth',3 );
        xlabel(handles.plot_pow,'Average Power');
    else
        WTamp = handles.amp_WT{signal_selected,1};
        handles.peak_value = max(WTamp(:))+.1;
        pcolor(handles.plot3d, handles.time_axis_us , handles.freqarr, WTamp);
        
        colorbar(handles.plot3d)
        plot(handles.plot_pow ,handles.amp_arr{signal_selected,1}, handles.freqarr,'-k','LineWidth',3 );
        xlabel(handles.plot_pow,'Average Amplitude');
    end
    
    colormap(handles.plot3d,handles.cmap);
    shading(handles.plot3d,'interp');
    
    if handles.calc_type == 1
        set(handles.plot3d,'yscale','log');
        set(handles.plot_pow,'yscale','log');
    end
    
    set(handles.plot3d,'ylim',[min(handles.freqarr) max(handles.freqarr)]);
    set(handles.plot3d,'xlim',[handles.time_axis_us(1) handles.time_axis_us(end)]);
    xlabel(handles.plot3d,'Time (s)');
    ylabel(handles.plot3d,'Frequency (Hz)');
    ylabel(handles.plot_pow,'Frequency (Hz)');
    ylim(handles.plot_pow,[min(handles.freqarr) max(handles.freqarr)]);
    set(handles.status,'String','Done Plotting');
end
grid(handles.plot3d,'on');
grid(handles.plot_pow,'off');
grid(handles.cum_avg,'off');

globalfontsize = 12; % Do not edit this line manually. See scripts/fontsize.py.

set(handles.plot3d,'FontUnits','points','FontSize',globalfontsize);
set(handles.plot_pow,'FontUnits','points','FontSize',globalfontsize);
set(handles.cum_avg,'FontUnits','points','FontSize',globalfontsize);

guidata(hObject,handles);

function intervals_Callback(hObject, eventdata, handles)
% Marks interval lines on plots
intervals = csv_to_mvar(get(handles.intervals,'String'));
intervals = sort(intervals);

% Clearing unmarked lines
child_handles = allchild(handles.wt_pane);

for i = 1:size(child_handles,1)
    if(strcmp(get(child_handles(i),'Type'),'axes'))
        axes_child = allchild(child_handles(i));
        for j = 1:size(axes_child,1)
            if strcmpi(get(axes_child(j),'Type'),'Line')
                line_style = get(axes_child(j),'linestyle');
                line_width = get(axes_child(j),'linewidth');
                if strcmp(line_style,'--') && line_width <= 1
                    delete(axes_child(j));
                end
            end
        end
        set(child_handles(i),'Ytickmode','auto','Xtickmode','auto');
        
    end
end

interval_selected = get(handles.signal_list,'Value');
hold(handles.cum_avg,'on');

if length(interval_selected)>1
    xl = get(handles.cum_avg,'ylim');
    for j = 1:size(intervals,2)
        x = [xl(1) xl(2)];
        z = ones(1,size(x,2));
        y = intervals(j)*ones(1,size(x,2));
        plot3(handles.cum_avg,y,x,z,'--k');
        xticks = get(handles.cum_avg,'xtick');
        xticks = unique(sort([xticks intervals]));
        set(handles.cum_avg,'xtick',xticks);
    end
    
elseif (interval_selected == size(handles.sig,1) + 1 ) %|| (~isfield(handles,'p'))
    xl = get(handles.cum_avg,'ylim');
    
    
    %
    for j = 1:size(intervals,2)
        x = [xl(1) xl(2)];
        z = ones(1,size(x,2));
        y = intervals(j)*ones(1,size(x,2));
        plot3(handles.cum_avg,y,x,z,'--k','HandleVisibility','off');
        xticks = get(handles.cum_avg,'xtick');
        xticks = unique(sort([xticks intervals]));
        set(handles.cum_avg,'xtick',xticks);
    end
    
    
else
    if interval_selected ~= size(handles.sig,1) + 1
        zval = 1;
        for i = 1:size(child_handles,1)
            if(strcmp(get(child_handles(i),'Type'),'axes') && strcmp(get(child_handles(i),'Visible'),'on'))
                hold(child_handles(i),'on');
                warning('off');
                xl = get(child_handles(i),'xlim');
                for j = 1:size(intervals,2)
                    
                    x = [xl(1) xl(2)];
                    z = ones(1,size(x,2));
                    z = z.*zval;
                    y = intervals(j)*ones(1,size(x,2));
                    plot3(child_handles(i),x,y,z,'--k');
                end
                yticks = get(child_handles(i),'ytick');
                yticks = unique(sort([yticks intervals]));
                set(child_handles(i),'Ytick',yticks);
                warning('on');
                hold(child_handles(i),'off');
            else
            end
            
            % Fix bug where xl has not been defined.
            if exist("xl", "var")
                h = child_handles(i);
                
                % Fix bug where limits of colorbar are set incorrectly.
                if ~isa(h, "matlab.graphics.illustration.ColorBar") && ~isa(h, "matlab.graphics.shape.internal.AnnotationPane")
                    try
                        set(h,'xlim',xl);
                    catch exception
                        disp("Error, possibly introduced by changes in a new version of MATLAB:");
                        disp(exception);
                    end
                end
            end
        end
    end
end

set(handles.plot_pow,'Yticklabel',[]);

function signal_list_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'r'
        signal_selected = get(handles.signal_list,'Value');
        if length(signal_selected)>1
            return;
        end
        list = get(handles.signal_list,'String');
        str = cell2mat(inputdlg('Enter the new signal name'));
        list{signal_selected,1} = str;
        set(handles.signal_list,'String',list);
end

%---------------------------Limits-----------------------------
function xlim_Callback(hObject, eventdata, handles)
% When the values of xlim are changed the graphs are updated
xl = csv_to_mvar(get(handles.xlim,'String'));
xlim(handles.time_series,xl);
xlim(handles.plot_pp,xl);
t = xl(2) - xl(1);
set(handles.length,'String',t);

function ylim_Callback(hObject, eventdata, handles)
% When the values of ylim are changed the graphs are updated
yl = csv_to_mvar(get(handles.ylim,'String'));
ylim(handles.time_series,yl);

%---------------------------Updating Value of limits Limits-----------------------------
function refresh_limits_Callback(hObject, eventdata, handles)
% Calculates limits of the plot
x = get(handles.time_series,'xlim');
xlim(handles.plot_pp,x);
t = x(2) - x(1);
x = strcat([num2str(x(1)),', ',num2str(x(2))]);

y = get(handles.time_series,'ylim');
y = strcat([num2str(y(1)),', ',num2str(y(2))]);

set(handles.xlim,'String',x);
set(handles.ylim,'String',y);
set(handles.length,'String',t);

% Update values in preprocess plot to reflect the
% preprocessed version of the x-limited signal.
preprocess_Callback(hObject, eventdata, handles);

% ---------------------------Zoom Updating--------------------------
function zoom_in_OffCallback(hObject, eventdata, handles)
% Refreshes the limit values right after the tool is deselected
x = get(handles.time_series,'xlim');
xlim(handles.plot_pp,x);
t = x(2) - x(1);
x = strcat([num2str(x(1)),', ',num2str(x(2))]);

y = get(handles.time_series,'ylim');
y = strcat([num2str(y(1)),', ',num2str(y(2))]);

set(handles.xlim,'String',x);
set(handles.ylim,'String',y);
set(handles.length,'String',t);

% -----------------------------Zoom Updating--------------------------
function zoom_out_OffCallback(hObject, eventdata, handles)
% Refreshes the limit values right after the tool is deselected
x = get(handles.time_series,'xlim');
xlim(handles.plot_pp,x);
t = x(2) - x(1);
x = strcat([num2str(x(1)),', ',num2str(x(2))]);

y = get(handles.time_series,'ylim');
y = strcat([num2str(y(1)),', ',num2str(y(2))]);

set(handles.xlim,'String',x);
set(handles.ylim,'String',y);
set(handles.length,'String',t);

function plot_type_SelectionChangeFcn(hObject, eventdata, handles)
% Deciding which plot type
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'power'
        handles.plot_type = 1;
    case 'amp'
        handles.plot_type = 2;
end

guidata(hObject,handles);
if ~isfield(handles,'p')
    xyplot_Callback(hObject, eventdata, handles)
elseif isempty(handles.p)
    xyplot_Callback(hObject, eventdata, handles)
else
    replot_Callback(hObject, eventdata, handles)
end
intervals_Callback(hObject, eventdata, handles)
guidata(hObject,handles);


function calc_type_SelectionChangedFcn(hObject, eventdata, handles)
% Deciding which type of calculation
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'wav'
        handles.calc_type = 1;
        list = {'Lognorm';'Morlet';'Bump';'';'';''};
        set(handles.wavelet_type,'String',list);
        set(handles.kaisera,'Enable','off')
    case 'four'
        handles.calc_type = 2;
        list = {'Gaussian';'Hann';'Blackman';'Exp';'Rect';'Kaiser'};
        set(handles.wavelet_type,'String',list);
        
end
drawnow;
guidata(hObject,handles);


function plot_TS_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.time_series, Fig);

globalfontsize = 12; % Do not edit this line manually. See scripts/fontsize.py.

set(ax,'Units', 'normalized', 'Position', [0.1,0.25,.85,.6],'FontUnits','points','FontSize',globalfontsize);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.3]);


function save_3dplot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.plot3d, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
colormap(ax,handles.cmap);
h = colorbar;
if handles.plot_type==1
    ylabel(h, 'Wavelet power')
else
    ylabel(h, 'Wavelet amplitude')
end


function save_avg_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.plot_pow, Fig);
view(90,-90);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);


function save_both_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.plot3d, Fig);
h = colorbar;
if handles.plot_type==1
    ylabel(h, 'Wavelet power')
else
    ylabel(h, 'Wavelet amplitude')
end
colormap(Fig,handles.cmap);
ax2 = copyobj(handles.plot_pow, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.07,0.2,.55,.7]);
set(ax2,'Units', 'normalized', 'Position', [0.75,0.2,.2,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
%ylabel(ax2,[])
set(Fig,'Units','normalized','Position', [0.2 0.2 0.6 0.5]);


function save_mm_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.cum_avg, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
if ~isfield(handles, "p") || isempty(handles.p)
    if ~isfield(handles, "leg1")
        handles.leg1={'Mean','Median'};
    end
    legend(ax,handles.leg1)
else
    legend(ax,handles.legstat)
end

function save_pow_arr_csv_Callback(hObject, eventdata, handles)
%Saves the avg power array in .csv format
[FileName,PathName] = uiputfile('.csv');
save_location = strcat(PathName,FileName);
pow_arr = cell2mat(handles.pow_arr);
csvwrite(save_location,pow_arr);

function save_amp_arr_csv_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('.csv');
save_location = strcat(PathName,FileName);
amp_arr = cell2mat(handles.amp_arr);
csvwrite(save_location,amp_arr);

function save_freqarr_csv_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('.csv');
save_location = strcat(PathName,FileName);
csvwrite(save_location,handles.freqarr');

function save_sig_pp_csv_Callback(hObject, eventdata, handles)
% Saves the preprocessed signal in .csv format
[FileName,PathName] = uiputfile('.csv');
save_location = strcat(PathName,FileName)
sig_pp = cell2mat(handles.sig_pp);
fs = handles.sampling_freq;
xl = get(handles.xlim,'String');
xl = csv_to_mvar(xl);
xl = xl.*fs;
xl(2) = min(xl(2),size(handles.sig,2));
xl(1) = max(xl(1),1);
sig_pp = sig_pp(:,xl(1):xl(2));
csvwrite(save_location,sig_pp);

function save_pow_arr_mat_Callback(hObject, eventdata, handles)
% Saves the avg power array in .mat format
[FileName,PathName] = uiputfile('.mat','Save Power Array as');
save_location = strcat(PathName,FileName)
powStruct.pow_arr = handles.pow_arr;
powStruct.freqarr = handles.freqarr';
save(save_location,'powStruct');

function save_amp_arr_mat_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('.mat','Save Amplitude Array as');
save_location = strcat(PathName,FileName)
ampStruct.amp_arr = handles.amp_arr;
ampStruct.freqarr = handles.freqarr';
save(save_location,'ampStruct');

function save_sig_pp_mat_Callback(hObject, eventdata, handles)
% Saves the preprocessed signal in .mat format
[FileName,PathName] = uiputfile('.mat','Save Preprocessed Signal as');
save_location = strcat(PathName,FileName)
sig_pp = cell2mat(handles.sig_pp);
fs = handles.sampling_freq;
xl = get(handles.xlim,'String');
xl = csv_to_mvar(xl);
xl = xl.*fs;
xl(2) = min(xl(2),size(handles.sig,2));
xl(1) = max(xl(1),1);
sig_pp = sig_pp(:,xl(1):xl(2));
save(save_location,'sig_pp');

function save_cut_ts_csv_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('.csv');
save_location = strcat(PathName,FileName);
sig = handles.sig;
fs = handles.sampling_freq;
xl = get(handles.xlim,'String');
xl = csv_to_mvar(xl);
xl = xl.*fs;
xl(2) = floor(min(xl(2),size(handles.sig,2)));
xl(1) = floor(max(xl(1),1));
sig = sig(:,xl(1):xl(2));
csvwrite(save_location,sig);

function save_cut_ts_mat_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('.mat','Save Cut Signal as');
save_location = strcat(PathName,FileName)
sig = handles.sig;
fs = handles.sampling_freq;
xl = get(handles.xlim,'String');
xl = csv_to_mvar(xl);
xl = xl.*fs;
xl(2) = min(xl(2),size(handles.sig,2));
xl(1) = max(xl(1),1);
sig = sig(:,xl(1):xl(2));
save(save_location,'sig');

function save_time_axis_mat_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('.mat','Save Cut Signal as');
save_location = strcat(PathName,FileName)
time_axis = get(handles.plot3d,'xlim');
fs = handles.sampling_freq;
time_axis = time_axis(1):1/fs:time_axis(2);
time_axis = time_axis(2:end);
save(save_location,'time_axis');

function save_time_axis_csv_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('.mat','Save Cut Signal as');
if FileName==0
    return;
else
end
save_location = strcat(PathName,FileName)
time_axis = get(handles.plot3d,'xlim');
fs = handles.sampling_freq;
time_axis = time_axis(1):1/fs:time_axis(2);
time_axis = time_axis(2:end);
csvwrite(save_location,time_axis);

function csv_save_Callback(hObject, eventdata, handles)
try
    [FileName,PathName] = uiputfile;
    save_location = strcat(PathName,FileName);
    save_location=save_location(1:end-3);
    xl = csv_to_mvar(get(handles.xlim,'String'));
    L=xl(2)*handles.wopt.fs - xl(1)*handles.wopt.fs;
    
    if strcmp(handles.wopt.Preprocess,'on')
        TFR_data.Preprocessed_data=handles.sig_pp;
    else
    end
    if handles.calc_type==1
        TFR_data.Analysis_type='Wavelet';
        TFR_data.Wavelet_type=handles.wopt.Wavelet;
    else
        TFR_data.Analysis_type='Windowed Fourier';
        TFR_data.Window_type=handles.wopt.Window;
    end
    if handles.plot_type==1
        avg_wt = cell2mat(handles.pow_arr);
        TFR_data.Power=avg_wt';
    else
        avg_wt = cell2mat(handles.amp_arr);
        TFR_data.Amplitude=avg_wt';
    end
    TFR_data.Frequency=handles.freqarr;
    TFR_data.Time=linspace(xl(1),xl(2),L);
    TFR_data.Sampling_frequency=handles.wopt.fs;
    TFR_data.fmax=handles.wopt.fmax;
    TFR_data.fmin=handles.wopt.fmin;
    TFR_data.fr=str2double(get(handles.central_freq,'String'));
    TFR_data.Preprocessing=handles.wopt.Preprocess;
    TFR_data.Cut_Edges=handles.wopt.CutEdges;
    
    if  isfield(handles,'p')
        TFR_data.group1_index=handles.g1;
        TFR_data.group2_index=handles.g2;
        TFR_data.group1=handles.gr1;
        TFR_data.group2=handles.gr2;
        TFR_data.p_value=handles.p;
        TFR_data.testtype=handles.ttypeS;
        x=get(handles.testinput,'String');
        TFR_data.testinput=x{handles.testin};
    else
    end
    
    
    data=csvsaving(TFR_data,handles.plot_type,handles.calc_type,handles);
    cell2csv([save_location,'_TFdata.csv'],data,',');
    
    clear data2
    D=TFR_data;
    if isfield(D,'Preprocessed_data')
        t=D.Time;
        pp=D.Preprocessed_data;
        S=size(pp,1);
        L=length(t);
        
        data2{1,1}='Time';
        for n=1:L
            data2{n+1,1}=t(n);
        end
        
        for j=1:S
            data2{1,1+j}=['Signal ',num2str(j)];
            for k=1:L
                data2{k+1,1+j}=pp{j}(k);
            end
        end
        cell2csv([save_location,'_ppdata.csv'],data2,',');
    else
    end
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end

function mat_save_Callback(hObject, eventdata, handles)
try
    [FileName,PathName] = uiputfile('.mat','Save data as');
    if FileName==0
        return;
    else
    end
    save_location = strcat(PathName,FileName);
    
    xl = csv_to_mvar(get(handles.xlim,'String'));
    L=xl(2)*handles.wopt.fs - xl(1)*handles.wopt.fs;
    
    if strcmp(handles.wopt.Preprocess,'on')
        TFR_data.Preprocessed_data=handles.sig_pp;
    else
    end
    if handles.calc_type==1
        TFR_data.Analysis_type='Wavelet';
        TFR_data.Wavelet_type=handles.wopt.Wavelet;
    else
        TFR_data.Analysis_type='Windowed Fourier';
        TFR_data.Window_type=handles.wopt.Window;
    end
    if handles.plot_type==1
        avg_wt = cell2mat(handles.pow_arr);
        TFR_data.Power=avg_wt';
    else
        avg_wt = cell2mat(handles.amp_arr);
        TFR_data.Amplitude=avg_wt';
    end
    if ~isempty(handles.currsig)
        TFR_data.Selected_sig=handles.currsig;
    else
    end
    TFR_data.Frequency=handles.freqarr;
    TFR_data.Time=linspace(xl(1),xl(2),L);
    TFR_data.Sampling_frequency=handles.wopt.fs;
    TFR_data.fmax=handles.wopt.fmax;
    TFR_data.fmin=handles.wopt.fmin;
    TFR_data.fr=str2double(get(handles.central_freq,'String'));
    TFR_data.Preprocessing=handles.wopt.Preprocess;
    TFR_data.Cut_Edges=handles.wopt.CutEdges;
    if  isfield(handles,'p')
        TFR_data.group1_index=handles.g1;
        TFR_data.group2_index=handles.g2;
        TFR_data.group1=handles.gr1;
        TFR_data.group2=handles.gr2;
        TFR_data.p_value=handles.p;
        TFR_data.testtype=handles.ttypeS;
        x=get(handles.testinput,'String');
        TFR_data.testinput=x{handles.testin};
    else
    end
    
    save(save_location,'TFR_data');
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end


function data=csvsaving(D,hp,hc,handles)

L=length(D.Frequency);

if hp==1;
    N=size(D.Power,2);
else
    N=size(D.Amplitude,2);
end

data=cell(L+13,(N*2)+1);
dstart=20;
data{1,1}='MODA v1.0 - Time-Frequency Analysis';
data{2,1}=date;
data{3,1}=[];
data{4,1}='PARAMETERS';
if hc==1
    data{5,1}='Analysis type';
    data{5,2}='Wavelet';
    data{6,1}='Wavelet type';
    data{6,2}=D.Wavelet_type;
else
    data{5,1}='Analysis type';
    data{5,2}='Windowed Fourier';
    data{6,1}='Window type';
    data{6,2}=D.Window_type;
end
data{7,1}='Sampling frequency (Hz)';
data{7,2}=D.Sampling_frequency;
data{8,1}='Maximum frequency (Hz)';
data{8,2}=D.fmax;
data{9,1}='Minimum frequency (Hz)';
data{9,2}=D.fmin;
data{10,1}='Frequency resolution';
data{10,2}=D.fr;
data{11,1}='Preprocessing';
data{11,2}=D.Preprocessing;

data{12,1}='Cut Edges';
data{12,2}=D.Cut_Edges;
data{13,1}='Time start (s)';
data{13,2}=min(D.Time);
data{14,1}='Time end (s)';
data{14,2}=max(D.Time);

if  isfield(D,'group1_index')
    
    data{15,1}='Group 1 signals';
    data{15,2}=D.group1_index;
    data{16,1}='Group 2 signals';
    data{16,2}=D.group2_index;
    data{17,1}='Test type';
    data{17,2}=D.testtype;
    data{18,1}='Test input';
    data{18,2}=D.testinput;
else
end



data{dstart,1}='Frequency';
for l=1:L;
    data{l+dstart,1}=D.Frequency(l);
end

if hp==1
    for j=1:N
        if isempty(handles.currsig)
            data{dstart,j+1}=['Power ',num2str(j)];
        else
            data{dstart,j+1}=['Power ',num2str(handles.currsig)];
        end
        for k=1:L
            data{k+dstart,j+1}=D.Power(k,j);
        end
        
    end
    
else
    
    for j=1:N
        if isempty(handles.currsig)
            data{dstart,j+1}=['Amplitude ',num2str(j)];
        else
            data{dstart,j+1}=['Amplitude ',num2str(handles.currsig)];
        end
        for k=1:L
            data{k+dstart,j+1}=D.Amplitude(k,j);
        end
        
    end
    
    
end

if  isfield(D,'group1_index')
    data{dstart,N+3}='p value (Group 1 vs group 2)';
    for k=1:L
        data{k+dstart,N+3}=D.p_value(k);
    end
else
end


% --------------------------------------------------------------------


% --------------------------------------------------------------------
function save_WT_coeff_Callback(hObject, eventdata, handles)
try
    [FileName,PathName] = uiputfile('.mat','Save data as');
    if FileName==0
        return;
    else
    end
    save_location = strcat(PathName,FileName);
    signal_selected = get(handles.signal_list, 'Value');
    xl = csv_to_mvar(get(handles.xlim,'String'));
    L=xl(2)*handles.wopt.fs - xl(1)*handles.wopt.fs;
    if strcmp(handles.wopt.Preprocess,'on')
        TFR_data.Preprocessed_data=handles.sig_pp{signal_selected,1};
    else
    end
    if handles.calc_type==1
        TFR_data.Analysis_type='Wavelet';
        TFR_data.Wavelet_type=handles.wopt.Wavelet;
    else
        TFR_data.Analysis_type='Windowed Fourier';
        TFR_data.Window_type=handles.wopt.Window;
    end
    TFR_data.TFcoefficients=handles.WT;
    
    TFR_data.Frequency=handles.freqarr;
    TFR_data.Time=linspace(xl(1),xl(2),L);
    TFR_data.Sampling_frequency=handles.wopt.fs;
    TFR_data.fmax=handles.wopt.fmax;
    TFR_data.fmin=handles.wopt.fmin;
    TFR_data.fr=str2double(get(handles.central_freq,'String'));%handles.wopt.f0;
    TFR_data.Preprocessing=handles.wopt.Preprocess;
    TFR_data.Cut_Edges=handles.wopt.CutEdges;
    
    save(save_location,'TFR_data');
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end


% --- Executes on selection change in wavelet_type.
function wavelet_type_Callback(hObject, eventdata, handles)


items = get(handles.wavelet_type,'String');
index_selected = get(handles.wavelet_type,'Value');
wtype = items{index_selected};
if strcmp(wtype,'Kaiser')
    set(handles.kaisera,'Enable','on')
else
    set(handles.kaisera,'Enable','off')
end


% --------------------------------------------------------------------
function resetGUI_Callback(hObject, eventdata, handles)

TimeFrequencyAnalysis;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

MODAclose(hObject,handles)


% --------------------------------------------------------------------
function save_session_Callback(hObject, eventdata, handles)
% hObject    handle to save_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MODAsave(handles)


% --------------------------------------------------------------------
function load_session_Callback(hObject, eventdata, handles)

handles=MODAload;



% % --- Executes on selection change in testtype.
% function testtype_Callback(hObject, eventdata, handles)
% replot_Callback(hObject, eventdata, handles)


% --- Executes on button press in replot.
function replot_Callback(hObject, eventdata, handles)
if isfield(handles,'amp_arr') || isfield(handles,'pow_arr')
    set(handles.status,'String','Plotting...')
    
    
    
    set(handles.save_WT_coeff,'Enable','off')
    cla(handles.cum_avg,'reset')
    cla(handles.plot3d,'reset')
    cla(handles.plot_pow,'reset')
    set(handles.plot3d,'Visible','off')
    set(handles.plot_pow,'Visible','off')
    hold(handles.cum_avg,'on')
    
    avgt=get(handles.avgtype,'Value');
    handles.ttype=get(handles.testtype,'Value');
    tlist=cellstr(get(handles.testtype,'String'));
    handles.ttypeS=tlist{handles.ttype};
    a=str2double(get(handles.alpha,'String'));
    handles.testin=get(handles.testinput,'Value');
    
    if handles.testin==2
        avg_wt = cell2mat(handles.pow_arr);
    else
        avg_wt = cell2mat(handles.amp_arr);
    end
    handles.g1=str2num(get(handles.group1,'String'));
    handles.g2=str2num(get(handles.group2,'String'));
    
    
    if isempty(handles.g1) || isempty(handles.g2)
        errordlg('Please enter signal numbers in group field','Error')
        return;
    end
    if length(handles.g1)<2 || length(handles.g2)<2
        errordlg('Group size must be larger than one for statistical analysis','Error')
        return;
    end
    
    handles.gr1=avg_wt(handles.g1,:);
    handles.gr2=avg_wt(handles.g2,:);
    
    if handles.ttype==2 && size(handles.gr1,1)~=size(handles.gr2,1)
        errordlg('Groups must be the same size for a paired test','Error')
        return;
    else
    end
    handles.h = waitbar(0,'Calculating statistics, please wait...');
    x=ones(1,length(handles.freqarr)).*a;
    handles.p=[];
    for j=1:length(handles.freqarr)
        if handles.ttype==2
            handles.p(j)=signrank(handles.gr1(:,j),handles.gr2(:,j));
        else
            handles.p(j)=ranksum(handles.gr1(:,j),handles.gr2(:,j));
        end
        
    end
    
    col=[0.9 0.9 0.9];
    
    if avgt==1
        
        plot(handles.cum_avg,handles.freqarr,median(handles.gr1),'Linewidth',2,'color',handles.linecol(1,:));
        plot(handles.cum_avg,handles.freqarr,median(handles.gr2),'Linewidth',2,'color',handles.linecol(2,:));
        axes(handles.cum_avg);
        val=fillsig(handles.freqarr,median(handles.gr1),median(handles.gr2),handles.p,x,'less',col);
        if val==1
            handles.legstat={'Group 1 median','Group 2 median',['p < ',num2str(a)]};
        else
            handles.legstat={'Group 1 median','Group 2 median'};
        end
        legend(handles.cum_avg,handles.legstat)
        if handles.testin==2
            ylabel(handles.cum_avg,'Median power')
        else
            ylabel(handles.cum_avg,'Median amplitude')
        end
        
        if handles.calc_type==1
            set(handles.cum_avg,'xscale','log')
        else
        end
        
        
    else
        plot(handles.cum_avg,handles.freqarr,mean(handles.gr1),'Linewidth',2,'color',handles.linecol(1,:))
        plot(handles.cum_avg,handles.freqarr,mean(handles.gr2),'Linewidth',2,'color',handles.linecol(2,:))
        axes(handles.cum_avg);
        val=fillsig(handles.freqarr,mean(handles.gr1),mean(handles.gr2),handles.p,x,'less',col);
        if val==1
            handles.legstat={'Group 1 median','Group 2 median',['p < ',num2str(a)]};
        else
            handles.legstat={'Group 1 median','Group 2 median'};
        end
        legend(handles.cum_avg,handles.legstat)
        if handles.testin==2
            ylabel(handles.cum_avg,'Mean power')
        else
            ylabel(handles.cum_avg,'Mean amplitude')
        end
        if handles.calc_type==1
            set(handles.cum_avg,'xscale','log')
        else
        end
    end
    
    
    set(handles.cum_avg,'xlim',[handles.freqarr(1) handles.freqarr(end)])
    box(handles.cum_avg,'on')
    title(handles.cum_avg,'Statistical comparison')
    xlabel(handles.cum_avg,'Frequency (Hz)')
    set(handles.save_3dplot,'Enable','off')
    set(handles.save_both_plot,'Enable','off')
    set(handles.save_avg_plot,'Enable','off')
    set(handles.save_mm_plot,'Enable','on')
    intervals_Callback(hObject, eventdata, handles)
    delete(handles.h)
    set(handles.status,'String','Done Plotting')
    guidata(hObject,handles);
else
end




% % --- Executes on selection change in avgtype.
% function avgtype_Callback(hObject, eventdata, handles)
% replot_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------



function max_freq_Callback(hObject, eventdata, handles)
preprocess_Callback(hObject, eventdata, handles)




function min_freq_Callback(hObject, eventdata, handles)
preprocess_Callback(hObject, eventdata, handles)
