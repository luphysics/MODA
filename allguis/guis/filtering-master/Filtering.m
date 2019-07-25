%Version 1.01
%**************************************************************************
%*************************** Filtering GUI ************************
%**************************************************************************
%---------------------------Credits---------------------------------------
% Wavelet Transform: Dmytro Iatsenko
% Ridge extraction: Dmytro Iatsenko
%
%----------------------------Documentation--------------------------------
%Reads a 1-D signal in either .mat or .csv format and displays it. 
%User can select the part of the signal he wants to use, and calculate wavelet
%tranform of that part. 
%Plots the Amplitude/Power surf plot and the average power plot over time. 
%Also contains save options for the graphs and data from wavelet transform.



function varargout = Filtering(varargin)
% FILTERING MATLAB code for Filtering.fig
%      FILTERING, by itself, creates a new FILTERING or raises the existing
%      singleton*.
%
%      H = FILTERING returns the handle to a new FILTERING or the handle to
%      the existing singleton*.
%
%      FILTERING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILTERING.M with the given input arguments.
%
%      FILTERING('Property','Value',...) creates a new FILTERING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Filtering_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Filtering_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help Filtering

% ------------------------------------------------------------------------
% Copyright notice for use of 'ginputc', downloaded on 18/09/2017 from
% Matlab FEX:
% Copyright (c) 2016, The MathWorks, Inc. 
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% * Redistributions of source code must retain the above copyright 
% notice, this list of conditions and the following disclaimer. 
% * Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in 
% the documentation and/or other materials provided with the distribution. 
% * In all cases, the software is, and all modifications and derivatives 
% of the software shall be, licensed to you solely for use in conjunction 
% with MathWorks products and service offerings.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% ------------------------------------------------------------------------

% Last Modified by GUIDE v2.5 20-Mar-2018 16:51:55
%*************************************************************************%
%                BEGIN initialization code - DO NOT EDIT                  %
%                ----------------------------------------                 %
%*************************************************************************%
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Filtering_OpeningFcn, ...
                   'gui_OutputFcn',  @Filtering_OutputFcn, ...
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

function Filtering_OpeningFcn(hObject, eventdata, handles, varargin)
% Executes when GUI is opened
handles=MODAsettings(hObject, handles);
handles.c=0; 
% Disable plotting and saving before data are loaded
set(handles.plot_TS,'Enable','off')
set(handles.save_3dplot,'Enable','off')
set(handles.save_both_plot,'Enable','off')
set(handles.save_avg_plot,'Enable','off')
set(handles.save_mm_plot,'Enable','off')
set(handles.save_filtered_sig_plot,'Enable','off')
set(handles.save_ridge_plot,'Enable','off') 
set(handles.save_phase_plot,'Enable','off') 
set(handles.All_filt_plot,'Enable','off') 
set(handles.save_fourier,'Enable','off') 
set(handles.save_csv,'Enable','off') 
set(handles.save_mat,'Enable','off') 
set(handles.save_session,'Enable','off')
set(handles.mark_interval,'Enable','off')
set(handles.add_interval,'Enable','off')
set(handles.filter_signal,'Enable','off')
set(handles.ridgecalc,'Enable','off')
handles.etype=2;


drawnow;
handles.output = hObject;
guidata(hObject, handles);

function varargout = Filtering_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function file_read_Callback(hObject, eventdata, handles)
[handles,A]=MODAreadcheck(handles);
if A==1
%-------Clearing axes and removing old data----------
    child_handles = allchild(handles.wt_pane);
    for i = 1:size(child_handles,1)
        if strcmp(get(child_handles(i),'type'),'axes')
            cla(child_handles(i),'reset');
            set(child_handles(i),'visible','off');
        end
    end
    
    %Remove frequency intervals from list, and freq input boxes.
    newItems = get(handles.interval_list,'String'); 
    newItems(:) = []; 
    set(handles.interval_list, 'String', newItems, 'Value', 1);
    set(handles.display_type,'Enable','off');


    set(handles.freq_1, 'String', []);
    set(handles.freq_2, 'String', []);
    if isfield(handles, 'freqarr');handles = rmfield(handles, 'freqarr');else end
    if isfield(handles, 'sig');handles = rmfield(handles, 'sig');else end
    if isfield(handles, 'sig_cut');handles = rmfield(handles, 'sig_cut');else end
    if isfield(handles, 'f1');handles = rmfield(handles, 'f1');else end    
    if isfield(handles, 'f2');handles = rmfield(handles, 'f2');else end     
    if isfield(handles, 'extract_phase');handles = rmfield(handles, 'extract_phase');else end
    if isfield(handles, 'extract_amp');handles = rmfield(handles, 'extract_amp');else end
    if isfield(handles, 'time_axis');handles = rmfield(handles, 'time_axis');else end
    if isfield(handles, 'pow_arr');handles = rmfield(handles, 'pow_arr');else end
    if isfield(handles, 'amp_arr');handles = rmfield(handles, 'amp_arr');else end
    if isfield(handles, 'pow_WT');handles = rmfield(handles, 'pow_WT');else end
    if isfield(handles, 'amp_WT');handles = rmfield(handles, 'amp_WT');else end
    if isfield(handles, 'WT');handles = rmfield(handles, 'WT');else end
    if isfield(handles, 'bands');handles = rmfield(handles, 'bands');else end
    if isfield(handles, 'wopt');handles = rmfield(handles, 'wopt');else end
    if isfield(handles, 'time_axis_us');handles = rmfield(handles, 'time_axis_us');else end
    if isfield(handles, 'sig_pp');handles = rmfield(handles, 'sig_pp');else end
    if isfield(handles, 'sampling_freq');handles = rmfield(handles, 'sampling_freq');else end
    if isfield(handles, 'peak_value');handles = rmfield(handles, 'peak_value');else end
    
    % Load data
    [handles,sig]=MODAread(handles,0);
    if sig==0
    else
    % Create signal list
    list = cell(size(sig,1)+1,1);
    list{1,1} = 'Signal 1';
    for i = 2:size(sig,1)
        list{i,1} = sprintf('Signal %d',i);
    end
    list{size(sig,1)+1,1} = sprintf('Average Plot (All)');
    set(handles.signal_list,'String',list); 

    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
         
    detrend_signal_Callback(hObject, eventdata, handles);%plots the detrended curve
%   xlabel(handles.time_series,'Time (s)','FontUnits','points','FontSize',10);
    
    set(handles.status,'String','Data loaded. Proceed with transform.');
    set(handles.transform,'Enable','on')    
    set(handles.plot_TS,'Enable','on')
    end

else
    return;
end

function refresh_limits_Callback(hObject, eventdata, handles)
%Calculates limits of the plot    
    x = get(handles.time_series,'xlim');
    y = get(handles.time_series,'ylim');
    
    xlim(handles.plot_pp,x);
    t = x(2) - x(1);
    
    xindex=x.*handles.sampling_freq;
    
    x = strcat([num2str(x(1)),', ',num2str(x(2))]);    
    y = strcat([num2str(y(1)),', ',num2str(y(2))]);
    %-------------------------
    
    xindex(2) = min(xindex(2),size(handles.sig,2));
    xindex(1) = max(xindex(1),1);
    handles.sig_cut=handles.sig(:,xindex(1):xindex(2));
    handles.time_axis_cut = handles.time_axis(xindex(1):xindex(2));
    handles.xl=[handles.time_axis_cut(1) handles.time_axis_cut(end)];
    
    %-------------------------
     
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);
    guidata(hObject,handles);
    
function detrend_signal_Callback(hObject, eventdata, handles)
% Preprocesses signal
ppstat=get(handles.preprocess,'Value');
if ppstat==2
    set(handles.plot_pp,'visible','on')
    set(handles.text38,'visible','on')
    cla(handles.plot_pp,'reset');
    sig_select=get(handles.signal_list,'Value');
    fmax = str2double(get(handles.max_freq,'String'));
    fmin = str2double(get(handles.min_freq,'String'));
    L = size(handles.sig_cut,2);
    N=size(handles.sig_cut);
    handles.sig_pp=NaN(N);
    
        for j = 1:size(handles.sig_cut,1)        
            sig = handles.sig_cut(j,:); 
            
            %Detrending
            X=(1:length(sig))'/handles.sampling_freq; XM=ones(length(X),4); 

            for pn=1:3 
                CX=X.^pn; 
                XM(:,pn+1)=(CX-mean(CX))/std(CX); 
            end
            sig = sig(:);
            
            w=warning('off','all'); 
            new_signal=sig-XM*(pinv(XM)*sig); 
            warning(w);

            %Filtering
            fx=fft(new_signal,L); % Fourier transform of a signal

            Nq=ceil((L+1)/2); 
            ff=[(0:Nq-1),-fliplr(1:L-Nq)]*handles.sampling_freq/L; 
            ff=ff(:); % frequencies in Fourier transform

            fx(abs(ff)<=max([fmin,handles.sampling_freq/L]) | abs(ff)>=fmax)=0; % filter signal in a chosen frequency domain
            handles.sig_pp(j,:) = ifft(fx)';
            
            % Plotting
            
            plot(handles.plot_pp,handles.time_axis_cut,handles.sig_cut(sig_select,:),'color',handles.linecol(1,:));
            hold(handles.plot_pp,'on');
            plot(handles.plot_pp,handles.time_axis_cut, handles.sig_pp(sig_select,:),'color',handles.linecol(2,:));
            legend(handles.plot_pp,{'Original','Pre-Processed'},'FontSize',8,'Location','Best','units','points');
            xlim(handles.plot_pp,[handles.time_axis_cut(1) handles.time_axis_cut(end)]);
            xlabel(handles.plot_pp,{'Time (s)'}); 
            
            drawnow;

        end   
else
    cla(handles.plot_pp,'reset');
    set(handles.plot_pp,'visible','off')
    set(handles.text38,'visible','off')
end
guidata(hObject,handles);

function preprocess_Callback(hObject, eventdata, handles)
% Executes when changing preprocessing preference
sig_select=get(handles.signal_list,'Value');
if sig_select==size(handles.sig_cut,1)+1
   set(handles.signal_list,'Value',1);
   detrend_signal_Callback(hObject, eventdata, handles)
   display_type_Callback(hObject, eventdata, handles)
else
    detrend_signal_Callback(hObject, eventdata, handles)
end

function signal_list_Callback(hObject, eventdata, handles)
% Executes when selected signal is changed

sig_select = get(handles.signal_list, 'Value');

if sig_select ~= size(handles.sig,1)+1 && length(sig_select)==1
    
    plot(handles.time_series,handles.time_axis_cut,handles.sig_cut(sig_select,:),'color',handles.linecol(1,:));
    xlim(handles.time_series,handles.xl)
    xlabel(handles.time_series,'Time (s)')
    
    detrend_signal_Callback(hObject, eventdata, handles)
    
    if isfield(handles,'freqarr')
        display_type_Callback(hObject, eventdata, handles)
    else
    end
else
    display_type_Callback(hObject, eventdata, handles)
end



function plot_type_SelectionChangeFcn(hObject, eventdata, handles)
% Executes when changing plot type
    switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
        case 'power'
            handles.plot_type = 1;
        case 'amp'
            handles.plot_type = 2;
    end
    
    guidata(hObject,handles); 
    
disp_select = get(handles.display_type,'Value');

if disp_select > 1
    return;
end
display_type_Callback(hObject, eventdata, handles)

% --- Executes when selected object is changed in calc_type.
function calc_type_SelectionChangedFcn(hObject, eventdata, handles)
% Executes when changing type of transform
switch get(eventdata.NewValue,'Tag')
    case 'wav'
        handles.calc_type=1;
        list={'Lognorm';'Morlet';'Bump';'';'';''};
        set(handles.wind_type,'String',list)
    case 'four'
        handles.calc_type=2;
        list={'Hann';'Gaussian';'Blackman';'Exp';'Rect';'Kaiser'};
        set(handles.wind_type,'String',list)
end
        
drawnow;
guidata(hObject,handles)


function wind_type_Callback(hObject, eventdata, handles)
% --- Executes on selection change in wind_type. Enables alpha input
% only when Kaiser window is selected.
wtypes = get(handles.wind_type,'String');
wselect = get(handles.wind_type,'Value');
wtype = wtypes{wselect}; 

if strcmp(wtype,'Kaiser')
    set(handles.kaisera,'Enable','on')
else
    set(handles.kaisera,'Enable','off')
end



function transform_Callback(hObject, eventdata, handles)
% --- Executes on button press in transform, and calculates the
% time-frequency representation of the data (WT or WFT).

set(handles.transform,'Enable','off')
set(handles.filter_signal,'Enable','off')
set(handles.ridgecalc,'Enable','off')

try
    
    
    set(handles.status,'String','Calculating Transform...');
    fmax=str2double(get(handles.max_freq,'String'));
    fmin=str2double(get(handles.min_freq,'String'));
    f0=str2double(get(handles.central_freq,'String')); handles.fc=f0;
    
    %% Do not allow f0<=0.4 for bump wavelet
    A=f0<=0.4;
    wtypes=get(handles.wind_type,'String');
    wselect=get(handles.wind_type,'Value');
    wtype=wtypes{wselect};
    B=strcmp(wtype,'Bump');
    
    if (A+0)+(B+0)==2
        errordlg('The bump wavelet requires that f0 > 0.4. Please enter a higher value.','Parameter Error');
        set(handles.transform,'Enable','on')
        return;
    end
    
    %% Do not allow maximum frequency to be higher than the Nyquist frequency (fs/2)
    if fmax>handles.sampling_freq/2
          errordlg(['Maximum frequency cannot be higher than the Nyquist frequency. Please enter a value less than or equal to ',num2str(fs/2),' Hz.'],'Parameter Error');
          set(handles.transform,'Enable','on')
          return;
    end
    
    %% Forces user to input minimum frequency for WFT, and changes resolution parameter according to fr/fmin, where fr is the user input resolution    
    if handles.calc_type==2 && isnan(fmin)
          errordlg(['Minimum frequency must be specified for WFT'],'Parameter Error');
          set(handles.transform,'Enable','on')
          return;
    elseif handles.calc_type==2
        handles.fc=f0/fmin;
    end
    
    %% Get user input alpha value if using Kaiser window
    if strcmp(wtype,'Kaiser')
        a=str2double(get(handles.kaisera,'String'));
        wtype = ['kaiser-',num2str(a)];        
    else
    end
    
    %% Prevents frequencies being too low
    if handles.calc_type==1
        if fmin<=1/(length(handles.sig_cut)/handles.sampling_freq)
          errordlg(['WT minimum frequency too low. To automatically calculate for minimum possible frequency leave "Min Freq" field blank.'],'Parameter Error'); 
          set(handles.wt_single,'Enable','on')
          set(handles.transform,'Enable','on')
          return;
        end
    else
    end
    
    % Preprocessing input
    x=get(handles.preprocess,'String'); ind=get(handles.preprocess,'Value'); ppselect=x{ind}; 
    
    % Cut edges input
    x=get(handles.cutedges,'String'); ind=get(handles.cutedges,'Value'); cutselect=x{ind};
    
    %% Downsample plotting of 2D results for speed
    N=length(handles.time_axis_cut);
    screensize=max(get(groot,'Screensize'));
    if N>screensize
        ds=floor(N/screensize);
    else
        ds=1;
    end
    
    handles.time_axis_ds=handles.time_axis_cut(1:ds:end);
    
    
    %% Create waitbar
    handles.h = waitbar(0,'Calculating transform...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(handles.h,'canceling',0)
    guidata(hObject,handles);
    
    
    %% Calculate the transform
    for p=1:size(handles.sig_cut,1)
        if getappdata(handles.h,'canceling') % Checks if user has pressed cancel on waitbar
            break;
        end
        
        set(handles.status,'String', sprintf('Calculating Transform of Signal %d/%d',p,size(handles.sig_cut,1)));
        wtwrapper; % Calls main calculation function
        handles.amp_WT{p}=abs(WT(:,1:ds:end));
        handles.pow_WT{p}=abs(WT(:,1:ds:end)).^2;
        handles.amp_av{p}=nanmean(handles.amp_WT{p},2); % Time-average
        handles.pow_av{p}=nanmean(handles.pow_WT{p},2);
        
        waitbar(p/length(handles.sig_cut),handles.h); % Update waitbar
        
    end
        
    set(handles.display_type,'Enable','on');
    set(handles.display_type,'Value',1);
    
    guidata(hObject,handles);
    
    display_type_Callback(hObject, eventdata, handles);

    delete(handles.h);
    set(handles.transform,'Enable','on')
    set(handles.filter_signal,'Enable','on')
    set(handles.ridgecalc,'Enable','on')
    set(handles.save_3dplot,'Enable','on')
    set(handles.save_both_plot,'Enable','on')
    set(handles.save_avg_plot,'Enable','on')
    set(handles.mark_interval,'Enable','on')
    set(handles.add_interval,'Enable','on')
    
    set(handles.file_read,'Enable','off')
    
catch e
    errordlg(e.message,'Error');
    set(handles.transform,'Enable','on')
    delete(handles.h);
    rethrow(e)
end

function display_type_Callback(hObject, eventdata, handles)
% Selecting what to display
disp_select=get(handles.display_type,'Value');
int_select=get(handles.interval_list,'Value');
sig_select=get(handles.signal_list,'Value');
extype=handles.etype;

% If single signal and time-frequency display are selected
if disp_select==1 && sig_select~=size(handles.sig_cut,1)+1
    
    if size(int_select,2)>1
        set(handles.interval_list,'max',1,'value',1);
    else
        set(handles.interval_list,'max',1);
    end   
    
    set(handles.save_3dplot,'Enable','on')
    set(handles.save_both_plot,'Enable','on')
    set(handles.save_avg_plot,'Enable','on')
    set(handles.save_mm_plot,'Enable','off')
    set(handles.save_filtered_sig_plot,'Enable','off')
    set(handles.save_ridge_plot,'Enable','off')
    set(handles.save_phase_plot,'Enable','off')
    set(handles.All_filt_plot,'Enable','off')
    set(handles.save_fourier,'Enable','off')    
    set(handles.fourier_scale,'visible','off')
    uistack(handles.plot3d,'top')
    uistack(handles.plot_pow,'top')
    set(handles.plot3d,'visible','on')
    set(handles.plot_pow,'visible','on')
    set(handles.mark_interval,'Enable','on')
    set(handles.add_interval,'Enable','on')
    cla(handles.cum_avg,'reset');
    set(handles.cum_avg,'visible','off');
    linkaxes([handles.amp_axis handles.phase_axis handles.freq_axis handles.time_series],'off');
    cla(handles.amp_axis,'reset');
    cla(handles.freq_axis,'reset');
    cla(handles.phase_axis,'reset');
    set(handles.amp_axis,'visible','off');
    set(handles.freq_axis,'visible','off');
    set(handles.phase_axis,'visible','off');
    cla(handles.fourier_plot,'reset');
    set(handles.fourier_plot,'visible','off');
    
    if handles.plot_type == 1     
          WTpow = handles.pow_WT{sig_select};
          handles.peak_value = max(WTpow(:))+.1;
          pcolor(handles.plot3d, handles.time_axis_ds , handles.freqarr, WTpow);                
          plot(handles.plot_pow, handles.pow_av{sig_select}, handles.freqarr,'-k','LineWidth',3 ,'color',handles.linecol(1,:));     
          xlabel(handles.plot_pow,'Average Power');
    else
          
          WTamp = handles.amp_WT{sig_select};
          handles.peak_value = max(WTamp(:))+.1;
          pcolor(handles.plot3d, handles.time_axis_ds , handles.freqarr, WTamp);         
          plot(handles.plot_pow ,handles.amp_av{sig_select}, handles.freqarr,'-k','LineWidth',3 ,'color',handles.linecol(1,:));
          xlabel(handles.plot_pow,'Average Amplitude');
    end
    
        colormap(handles.plot3d,handles.cmap);
        shading(handles.plot3d,'interp');
        set(handles.plot3d,'yscale','log');
        set(handles.plot_pow,'yscale','log');
        
        if handles.calc_type==2
            set(handles.plot3d,'yscale','linear');
            set(handles.plot_pow,'yscale','linear');
        else
            
        end
        
        set(handles.plot3d,'ylim',[min(handles.freqarr) max(handles.freqarr)]);%making the axes tight
        set(handles.plot3d,'xlim',[handles.time_axis_ds(1) handles.time_axis_ds(end)]);%making the axes tight
        xlabel(handles.plot3d,'Time (s)');
        ylabel(handles.plot3d,'Frequency (Hz)');    
        ylabel(handles.plot_pow,'Frequency (Hz)');    
        ylim(handles.plot_pow,[min(handles.freqarr) max(handles.freqarr)]);
        set(handles.status,'String','Done Plotting');
        
% If average plot and time-frequency display are selected    
elseif disp_select==1 && sig_select==size(handles.sig_cut,1)+1
     if isfield(handles,'freqarr')
        set(handles.mark_interval,'Enable','on')
        set(handles.add_interval,'Enable','on')
        set(handles.save_3dplot,'Enable','off')
        set(handles.save_both_plot,'Enable','off')
        set(handles.save_avg_plot,'Enable','off')
        set(handles.save_mm_plot,'Enable','on')
        set(handles.save_filtered_sig_plot,'Enable','off')
        set(handles.save_ridge_plot,'Enable','off')
        set(handles.save_phase_plot,'Enable','off')
        set(handles.All_filt_plot,'Enable','off')
        set(handles.save_fourier,'Enable','off') 
        cla(handles.plot3d,'reset');
        cla(handles.plot_pow,'reset');
        cla(handles.cum_avg,'reset');
        set(handles.plot3d,'visible','off');
        set(handles.plot_pow,'visible','off');   
        set(handles.cum_avg,'visible','on');
        hold(handles.cum_avg,'on');
        uistack(handles.cum_avg,'top');
    
       
        if(handles.plot_type == 1)      
            plot(handles.cum_avg, handles.freqarr, mean(cell2mat(handles.pow_av),2),'-','Linewidth',3,'color',handles.linecol(1,:));
            plot(handles.cum_avg, handles.freqarr, median(cell2mat(handles.pow_av),2),'--','Linewidth',3,'color',handles.linecol(2,:));
            ylabel(handles.cum_avg,'Average Power');
            xlabel(handles.cum_avg,'Frequency (Hz)');
        else
            plot(handles.cum_avg, handles.freqarr, mean(cell2mat(handles.amp_av),2),'-','Linewidth',3,'color',handles.linecol(1,:));
            plot(handles.cum_avg, handles.freqarr, median(cell2mat(handles.amp_av),2),'--','Linewidth',3,'color',handles.linecol(2,:));
            ylabel(handles.cum_avg,'Average Amplitude');
            xlabel(handles.cum_avg,'Frequency (Hz)');
            
        end
        
        set(handles.cum_avg,'xscale','log');
        
        if handles.calc_type==2
            set(handles.cum_avg,'xscale','linear');
        else
            
        end
        
        handles.leg1={'Mean','Median'};
        legend(handles.cum_avg,handles.leg1)
        xlim(handles.cum_avg,[handles.freqarr(1) handles.freqarr(end)])
        else
    end
 
% If 'Bands' and a single signal are selected
elseif disp_select==2
if ~isfield(handles,'bands') && ~isfield(handles,'recon')
    return;
else
    if sig_select==size(handles.sig_cut,1)+1
        set(handles.signal_list,'Value',1)
        sig_select=1;
    else
    end
    
    if isempty(int_select) 
            return;
    end
        
    set(handles.save_3dplot,'Enable','off')
    set(handles.save_both_plot,'Enable','off')
    set(handles.save_avg_plot,'Enable','off')
    set(handles.save_mm_plot,'Enable','off')
    set(handles.save_filtered_sig_plot,'Enable','on')
    set(handles.save_ridge_plot,'Enable','on')
    set(handles.save_phase_plot,'Enable','on')
    set(handles.All_filt_plot,'Enable','on')
    set(handles.save_fourier,'Enable','off')  
    set(handles.fourier_scale,'visible','off');
    list = get(handles.interval_list,'String');
    
    set(handles.interval_list,'max',size(list,1));
    cla(handles.plot3d,'reset');
    cla(handles.plot_pow,'reset');
    cla(handles.cum_avg,'reset');
    cla(handles.fourier_plot,'reset');
    set(handles.plot3d,'visible','off');
    set(handles.plot_pow,'visible','off');
    set(handles.cum_avg,'visible','off');
    set(handles.fourier_plot,'visible','off');
 
    uistack(handles.amp_axis,'top');
    uistack(handles.phase_axis,'top');
    uistack(handles.freq_axis,'top');
    cla(handles.amp_axis,'reset');
    cla(handles.phase_axis,'reset');
    cla(handles.freq_axis,'reset');
    set(handles.amp_axis,'visible','on');
    set(handles.phase_axis,'visible','on');
    set(handles.freq_axis,'visible','on');
    if ~isfield(handles,'bands') && ~isfield(handles,'bands_iphi')
            return;
    end
        
    
    hold(handles.amp_axis,'on');
    hold(handles.phase_axis,'on');
    hold(handles.freq_axis,'on');
    
    fb=get(handles.interval_list,'String');
    fb=cell2mat(fb);

    N=size(fb,1);

    for j=1:N
        handles.f1{j}=fb(j,1:4);
        handles.f2{j}=fb(j,10:13);
    end
    
    if extype==1
        for j=1:size(int_select,2)
            plot(handles.amp_axis,handles.time_axis_cut,handles.recon{sig_select,int_select(j)},'color',handles.linecol(j,:),'linewidth',handles.line2width);
            handles.leg3{j}=[handles.f1{int_select(j)},' - ',handles.f2{int_select(j)},' Hz'];
            legend(handles.amp_axis,handles.leg3,'FontSize',10,'Orientation','Vertical','position',[0.89 0.5 0.05 0.05])
            plot(handles.phase_axis, handles.time_axis_cut,handles.bands_iphi{sig_select,int_select(j)},'color',handles.linecol(j,:),'linewidth',handles.line2width);
            
        
        
        if size(int_select,2)==1
            if(handles.plot_type == 1)      
                WTpow = handles.pow_WT{sig_select};
                handles.peak_value = max(WTpow(:))+.1;
                pcolor(handles.freq_axis, handles.time_axis_ds , handles.freqarr, WTpow(1:end,1:end));                     
            else 
                WTamp = handles.amp_WT{sig_select};
                handles.peak_value = max(WTamp(:))+.1;
                pcolor(handles.freq_axis, handles.time_axis_ds , handles.freqarr, WTamp(1:end,1:end));         
            end
            
            colormap(handles.cmap);
            shading(handles.freq_axis,'interp');
            ylim(handles.freq_axis,csv_to_mvar(list{int_select,1}));
             plot(handles.freq_axis, handles.time_axis_cut,handles.bands_freq{sig_select,int_select},'color',handles.linecol(j,:),'linewidth',handles.line2width);
        else
            
            shading(handles.freq_axis,'interp');
            plot(handles.freq_axis, handles.time_axis_cut,handles.bands_freq{sig_select,int_select(j)},'color',handles.linecol(j,:),'linewidth',handles.line2width);
        end
        end
        
        guidata(hObject,handles);
        if handles.calc_type==2
                set(handles.freq_axis,'yscale','linear');
        else
                set(handles.freq_axis,'yscale','log');
        end
        linkaxes([handles.amp_axis handles.phase_axis handles.freq_axis handles.time_series],'x');
        xlim(handles.amp_axis,handles.xl);
        xlim(handles.phase_axis,handles.xl);       
        xlabel(handles.phase_axis,'Time (s)');
        ylabel(handles.phase_axis,'Phase');
        ylabel(handles.amp_axis,'Filtered Signal');
        ylabel(handles.freq_axis,'Frequency (Hz)');
        set(handles.phase_axis,'yticklabel',{'-\pi','-0.5\pi','0', '0.5\pi', '\pi'},'ytick',[-pi, -0.5*pi, 0, 0.5*pi, pi],'fontunits','normalized');
        set(handles.amp_axis,'fontunits','points','fontsize',10);
        
        
    elseif extype==2
        for j=1:size(int_select,2)
            plot(handles.amp_axis,handles.time_axis_cut,handles.bands{sig_select,int_select(j)},'color',handles.linecol(j,:),'linewidth',handles.line2width);
            handles.leg3{j}=[handles.f1{int_select(j)},' - ',handles.f2{int_select(j)},' Hz']; 
            legend(handles.amp_axis,handles.leg3,'FontSize',10,'Orientation','Vertical','position',[0.89 0.5 0.05 0.05])
            plot(handles.phase_axis, handles.time_axis_cut, handles.extract_phase{sig_select,int_select(j)},'color',handles.linecol(j,:),'linewidth',handles.line2width);
            plot(handles.freq_axis,handles.time_axis_cut,handles.extract_amp{sig_select,int_select(j)},'color',handles.linecol(j,:),'linewidth',handles.line2width);
        end
        guidata(hObject,handles);
        linkaxes([handles.amp_axis handles.phase_axis handles.freq_axis handles.time_series],'x');
        xlim(handles.amp_axis,handles.xl);
        xlim(handles.phase_axis,handles.xl);       
        xlabel(handles.phase_axis,'Time (s)');
        ylabel(handles.phase_axis,'Phase');
        ylabel(handles.amp_axis,'Filtered Signal');
        ylabel(handles.freq_axis,'Amplitude');
        set(handles.amp_axis,'fontunits','points','fontsize',10);
    else
    end
end
% If Fourier plot is selected
elseif disp_select==3
    if ~isfield(handles,'bands') && ~isfield(handles,'recon')
    return;
    else
    set(handles.save_3dplot,'Enable','off')
    set(handles.save_both_plot,'Enable','off')
    set(handles.save_avg_plot,'Enable','off')
    set(handles.save_mm_plot,'Enable','off')
    set(handles.save_filtered_sig_plot,'Enable','off')
    set(handles.save_ridge_plot,'Enable','off')
    set(handles.save_phase_plot,'Enable','off')
    set(handles.All_filt_plot,'Enable','off')
    set(handles.save_fourier,'Enable','on') 
    set(handles.fourier_scale,'visible','on');
    list = get(handles.interval_list,'String');
    set(handles.interval_list,'max',size(list,1));
    linkaxes([handles.amp_axis handles.phase_axis handles.freq_axis handles.time_series],'off');    
    cla(handles.plot3d,'reset');
    cla(handles.plot_pow,'reset');
    cla(handles.cum_avg,'reset');
    cla(handles.amp_axis,'reset');
    cla(handles.fourier_plot,'reset');
    set(handles.plot3d,'visible','off');
    set(handles.plot_pow,'visible','off');
    set(handles.amp_axis,'visible','off');
    set(handles.freq_axis,'visible','off');
    set(handles.phase_axis,'visible','off');   
    set(handles.cum_avg,'visible','off');
    set(handles.fourier_plot,'visible','on');
    hold(handles.fourier_plot,'on');
    uistack(handles.fourier_plot,'top');
    hold(handles.fourier_plot,'on');

    % Preprocessing input
    x=get(handles.preprocess,'String'); ind=get(handles.preprocess,'Value'); ppselect=x{ind};
    
    if strcmp(ppselect,'on')
        [ftorig, ft_freq] = Fourier(handles.sig_pp(sig_select,:),handles.sampling_freq);
    else
        [ftorig, ft_freq] = Fourier(handles.sig(sig_select,:),handles.sampling_freq);
    end
    
    plot(handles.fourier_plot,ft_freq,ftorig,'linewidth',2,'color',handles.linecol(1,:));
    set(handles.fourier_plot,'fontunits','points','fontsize',10,'visible','on','xscale','log','yscale','log','xlim',[handles.freqarr(1) handles.freqarr(end)]);
    
    extype = handles.etype;
    handles.leg2={'Original'}; 
    fb=get(handles.interval_list,'String');
    fb=cell2mat(fb);

    N=size(fb,1);

    for j=1:N
        handles.f1{j}=fb(j,1:4);
        handles.f2{j}=fb(j,10:13);
    end
    
    legend(handles.fourier_plot,handles.leg2,'FontSize',10)
    if isfield(handles,'bands') && extype==2
        for j = 1:size(int_select,2)        
            [ft, ft_freq] = Fourier(handles.bands{sig_select,int_select(j)},handles.sampling_freq);            
            plot(handles.fourier_plot,ft_freq,ft,'Linewidth',2,'color',handles.linecol(j+1,:));            
            handles.leg2{j+1}=[handles.f1{int_select(j)},' - ',handles.f2{int_select(j)},' Hz']; 
            legend(handles.fourier_plot,handles.leg2,'FontSize',10)
        end
        
        
        
        yl = get(handles.fourier_plot,'ylim'); 
%         for j = 1:size(int_select,2)   
%             fl = csv_to_mvar(list{int_select(j),1});            
%             x = fl(1)*[1 1];
%             plot(handles.fourier_plot,x,yl,'-k');
%             x = fl(2)*[1 1];
%             plot(handles.fourier_plot,x,yl,'-k');
%         end
    end     
    
    if isfield(handles,'bands_iphi') && extype==1
        for j = 1:size(int_select,2)        
            [ft, ft_freq] = Fourier(handles.recon{sig_select,int_select(j)},handles.sampling_freq);            
            plot(handles.fourier_plot,ft_freq,ft,'Linewidth',2,'color',handles.linecol(j+1,:));            
            handles.leg2{j+1}=[handles.f1{int_select(j)},' - ',handles.f2{int_select(j)},' Hz']; 
            legend(handles.fourier_plot,handles.leg2,'FontSize',10)
            set(handles.fourier_plot,'ylim',[min(abs(ft)) max(abs(ft))]);
        end    
    end
    
    %fourier_scale_Callback(hObject, eventdata, handles)
    xlabel(handles.fourier_plot,'Frequency (Hz)','fontunits','points','fontsize',10)
    ylabel(handles.fourier_plot,'FT Power','fontunits','points','fontsize',10)
    
    guidata(hObject,handles);
    
    end
    
    
    x=get(handles.fourier_plot,'ylim');
    set(handles.fourier_plot,'ylim',[x(1) max(ftorig)]);
    
    x= get(handles.fourier_scale,'Value');
    if x==1
        set(handles.fourier_plot,'xscale','log','yscale','log')
    else
        set(handles.fourier_plot,'xscale','linear','yscale','linear')
    end
    
    for j = 1:size(int_select,2)   
        yl = get(handles.fourier_plot,'ylim'); 
        fl = csv_to_mvar(list{int_select(j),1});            
        x = fl(1)*[1 1];
        plot(handles.fourier_plot,x,yl,'-k');
        x = fl(2)*[1 1];
        plot(handles.fourier_plot,x,yl,'-k');
    end
    
end

function max_freq_Callback(hObject, eventdata, handles)
    preprocess_Callback(hObject, eventdata, handles)
    
function min_freq_Callback(hObject, eventdata, handles)
    preprocess_Callback(hObject, eventdata, handles)
    
    % --- Executes on button press in ridgecalc.
function ridgecalc_Callback(hObject, eventdata, handles)
handles.etype=1;
handles=MODAridge_filter(hObject,eventdata,handles);
guidata(hObject,handles);

display_type_Callback(hObject, eventdata, handles)

function filter_signal_Callback(hObject, eventdata, handles)
% Executes when user presses 'Filter Signal'
handles.etype=2;
handles=MODAridge_filter(hObject,eventdata,handles);
guidata(hObject,handles);
display_type_Callback(hObject, eventdata, handles)

    
function mark_interval_Callback(hObject, eventdata, handles)
% Executes when 'Mark Interval' is pressed
disp_select = get(handles.display_type,'Value');
sig_select=get(handles.signal_list,'Value');

if disp_select ~= 1
    return;
end



if any(sig_select == size(handles.sig,1)+1)
    
    child_handles = allchild(handles.cum_avg);
        for i = 1:size(child_handles,1)   
            line_style = get(child_handles(i),'linestyle');
            line_size = get(child_handles(i),'linewidth');
            if(strcmp(get(child_handles(i),'Type'),'line') && strcmp(line_style,'--')) && str2double(line_size)<1
                    delete(child_handles(i))
            end
        end
        
        grid(handles.cum_avg,'on');
        hold(handles.cum_avg,'on');  
        
        [f,~] = ginput(1);
        set(handles.freq_1,'String',f);
        ylavg = get(handles.cum_avg,'ylim');
        line([f f],ylavg,'Color','k','LineStyle','--');
        
        [f,~] = ginput(1);
        set(handles.freq_2,'String',f);
        ylavg = get(handles.cum_avg,'ylim');
        line([f f],ylavg,'Color','k','LineStyle','--');
        
else
    xl3d = get(handles.plot3d,'xlim');
    xlpow = get(handles.plot_pow,'xlim');  
    z = [1 1];
    clear_axes_lines(handles.plot3d);
    child_handles = allchild(handles.plot_pow);
        for i = 1:size(child_handles,1)   
            line_style = get(child_handles(i),'linestyle');
            if(strcmp(get(child_handles(i),'Type'),'line') && strcmp(line_style,'--')) 
            delete(child_handles(i))
            end
        end
        
        grid(handles.plot_pow,'off');
        grid(handles.plot3d,'on');
        hold(handles.plot3d,'on');     
        hold(handles.plot_pow,'on');
        xlim(handles.plot_pow,xlpow);
        
        [~,f] = ginput(1);
        set(handles.freq_1,'String',f);
        plot3(handles.plot3d,xl3d,[f f],z,'--k');
        plot(handles.plot_pow,xlpow,[f f],'--k');
        
        [~,f] = ginput(1);
        set(handles.freq_2,'String',f);
        plot3(handles.plot3d,xl3d,[f f],z,'--k');
        plot(handles.plot_pow,xlpow,[f f],'--k');
        
end


hold(handles.plot3d,'off');  
hold(handles.plot_pow,'off');
hold(handles.cum_avg,'off');


function add_interval_Callback(hObject, eventdata, handles)
% Executes when 'Add interval' is pressed
handles.c=handles.c+1;
f1 = str2double(get(handles.freq_1,'String'));
f2 = str2double(get(handles.freq_2,'String'));

freqmin=min(handles.freqarr);
freqmax=max(handles.freqarr);

if f1<freqmin || f1>freqmax
      errordlg('Selected frequencies are outside the allowable range','Parameter Error');
      return;
end

if f2<freqmin || f2>freqmax
      errordlg('Selected frequencies are outside the allowable range','Parameter Error');
      return;
end    

fl = sprintf('%f,%f',min(f1,f2),max(f1,f2));
list = get(handles.interval_list,'String');
list{end+1,1} = fl;
set(handles.interval_list,'String',list);
guidata(hObject,handles);
drawnow;

function interval_list_Callback(hObject, eventdata, handles)
display_type_Callback(hObject, eventdata, handles)

function fourier_scale_Callback(hObject, eventdata, handles)
    
    display_type_Callback(hObject, eventdata, handles)
    

%% Plotting functions
function plot_TS_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.time_series, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.25,.85,.6],'FontUnits','points','FontSize',10);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.3]);


function save_3dplot_Callback(hObject, eventdata, handles)
%Saves the 3d plot
Fig = figure;
ax = copyobj(handles.plot3d, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
colormap(ax,handles.cmap);
colorbar

function save_avg_plot_Callback(hObject, eventdata, handles)
%Saves the power plot
Fig = figure;
ax = copyobj(handles.plot_pow, Fig);
view(90,-90);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);

function save_both_plot_Callback(hObject, eventdata, handles)
%Saves the 3D and power plot
Fig = figure;
ax1 = copyobj(handles.plot3d, Fig);
ax2 = copyobj(handles.plot_pow, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.2,.55,.7]);
set(ax2,'Units', 'normalized', 'Position', [0.7,0.2,.25,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
ylabel(ax2,[])
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
colorbar
colormap(Fig,handles.cmap);

function save_mm_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.cum_avg, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
legend(ax,handles.leg1,'FontSize',10)

function save_filtered_sig_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.amp_axis, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
ylabel(ax,'Filtered Signal')
xlabel(ax,'Time (s)')
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
legend(ax,handles.leg3,'Orientation','Horizontal')

function save_ridge_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.freq_axis, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
xlabel('Time (s)')
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);

colormap(ax,handles.cmap)

interval_selected = get(handles.interval_list,'Value');
    if size(interval_selected,2)>1
        legend(ax,handles.leg3,'orientation','horizontal')
    else
    end

function save_phase_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.phase_axis, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
% set(ax,'FontSize',0.05)
% ylabel(ax,'Phase','FontSize',0.05)
% xlabel(ax,'Time (s)','FontSize',0.05)
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
legend(ax,handles.leg3,'orientation','horizontal')



function All_filt_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.amp_axis, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.7,.85,.25]);
ylabel(ax,'Filtered Signal')
legend(ax,handles.leg3,'orientation','horizontal')

ax2 = copyobj(handles.freq_axis, Fig);
set(ax2,'Units', 'normalized', 'Position', [0.1,0.4,.85,.25]);

if handles.etype==2
    ylabel(ax2,'Amplitude')
else
    ylabel(ax2,'Frequency (Hz)')
    colormap(handles.cmap)
end

ax3 = copyobj(handles.phase_axis, Fig);
set(ax3,'Units', 'normalized', 'Position', [0.1,0.1,.85,.25]);
ylabel(ax3,'Phase')
xlabel(ax3,'Time (s)')
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
colormap(handles.cmap)

function save_fourier_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.fourier_plot, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
legend(ax,handles.leg2)

%% Saving
function save_csv_Callback(hObject, eventdata, handles)
try
% [FileName,PathName] = uiputfile('.csv','Save as');
% save_location = strcat(PathName,FileName)

curr=pwd;

[FileName,PathName] = uiputfile('.csv','Save as');
if FileName==0
    return;
else
end
cd(PathName)

foldername=[FileName(1:end-4)];
mkdir(foldername)


xl = csv_to_mvar(get(handles.xlim,'String'));
L=xl(2)*handles.wopt.fs - xl(1)*handles.wopt.fs;
list = get(handles.interval_list,'String');

Filtered_data.Sampling_frequency = handles.wopt.fs;
Filtered_data.Time=linspace(xl(1),xl(2),L);
Filtered_data.Freq_bands=list;
if handles.etype==1
    Filtered_data.Filter_type='Ridge Extraction';
else
    Filtered_data.Filter_type='Butterworth Filter';
end

if handles.etype==1
    
if handles.calc_type==1
    Filtered_data.Analysis_type='Wavelet';
    Filtered_data.Wavelet_type=handles.wopt.Wavelet;
else
    Filtered_data.Analysis_type='Windowed Fourier';
    Filtered_data.Window_type=handles.wopt.Window;
end
    Filtered_data.Preprocessing=handles.wopt.Preprocess;
    %Filtered_data.Cut_Edges=handles.wopt.CutEdges;
    Filtered_data.Frequency_resolution=str2double(get(handles.central_freq,'String'));
    Filtered_data.Ridge_recon=handles.recon;
    Filtered_data.Ridge_frequency=handles.bands_freq;
    Filtered_data.Ridge_amplitude=handles.bands_iamp;
    Filtered_data.Ridge_phase=handles.bands_iphi;
else
    Filtered_data.Filtered_sigs=handles.bands;
    Filtered_data.Filtered_phases=handles.extract_phase;
    Filtered_data.Filtered_amplitudes=handles.extract_amp;
end

data1=MODAcsvsave(Filtered_data,1);
data2=MODAcsvsave(Filtered_data,2);
data3=MODAcsvsave(Filtered_data,3);

if handles.etype==1
    cell2csv([foldername,'\extracted_modes.csv'],data1,',');
    cell2csv([foldername,'\ridge_frequencies.csv'],data2,',');
    cell2csv([foldername,'\ridge_amplitudes.csv'],data3,',');
else    
    cell2csv([foldername,'\filtered_signals.csv'],data1,',');
    cell2csv([foldername,'\filtered_phases.csv'],data2,',');
    cell2csv([foldername,'\filtered_amplitudes.csv'],data3,',');
end
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end


function save_mat_Callback(hObject, eventdata, handles)
try
[FileName,PathName] = uiputfile('.mat','Save as');
if FileName==0
    return;
else
end
save_location = strcat(PathName,FileName)

xl = csv_to_mvar(get(handles.xlim,'String'));
L=xl(2)*handles.wopt.fs - xl(1)*handles.wopt.fs;
list = get(handles.interval_list,'String');

Filtered_data.Sampling_frequency = handles.wopt.fs;
Filtered_data.Time=linspace(xl(1),xl(2),L);

Filtered_data.Freq_bands=list;
if handles.etype==1
    Filtered_data.Filter_type='Ridge Extraction';
else
    Filtered_data.Filter_type='Butterworth Filter';
end

if handles.etype==1
    
if handles.calc_type==1
    Filtered_data.Analysis_type='Wavelet';
    Filtered_data.Wavelet_type=handles.wopt.Wavelet;
else
    Filtered_data.Analysis_type='Windowed Fourier';
    Filtered_data.Window_type=handles.wopt.Window;
end
    Filtered_data.Preprocessing=handles.wopt.Preprocess;
    %Filtered_data.Cut_Edges=handles.wopt.CutEdges;
    Filtered_data.Frequency_resolution=str2double(get(handles.central_freq,'String'));%handles.wopt.f0;
    Filtered_data.Ridge_recon=handles.recon;
    Filtered_data.Ridge_frequency=handles.bands_freq;
    Filtered_data.Ridge_amplitude=handles.bands_iamp;
    Filtered_data.Ridge_phase=handles.bands_iphi;
else
    Filtered_data.Filtered_sigs=handles.bands;
    Filtered_data.Filtered_phases=handles.extract_phase;
    Filtered_data.Filtered_amplitudes=handles.extract_amp;
end
save(save_location,'Filtered_data');
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

function resetGUI_Callback(hObject, eventdata, handles)
% Executes when requested new workspace
Filtering;


function figure1_CloseRequestFcn(hObject, eventdata, handles)
% --- Executes when user attempts to close figure1.
MODAclose(hObject,handles)

function save_session_Callback(hObject, eventdata, handles)

MODAsave(handles)


% --------------------------------------------------------------------
function load_session_Callback(hObject, eventdata, handles)

handles=MODAload;


function interval_list_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'delete'
        
        interval_selected = get(handles.interval_list,'Value');
        if min(interval_selected)>1
            set(handles.interval_list,'Value',min(interval_selected)-1);
        else
            set(handles.interval_list,'Value',1);
        end
        list = get(handles.interval_list,'String');
        list(interval_selected,:) = [];
        set(handles.interval_list,'String',list);   
        n=1:handles.c;
        ne=n(1:end ~=interval_selected);
        
        if isfield(handles,'bands')
            handles.bands=handles.bands(:,ne);
            handles.extract_phase=handles.extract_phase(:,ne);
            handles.extract_amp=handles.extract_amp(:,ne);
        end
        
        if isfield(handles,'recon')
            handles.bands_iamp=handles.bands_iamp(:,ne);
            handles.bands_iphi=handles.bands_iphi(:,ne);
            handles.bands_freq=handles.bands_freq(:,ne);
            handles.recon=handles.recon(:,ne);
        end
        
        
%         handles.bands(:,interval_selected) = [];
handles.c=handles.c-1;
%         guidata(hObject,handles);
         interval_list_Callback(hObject, eventdata, handles)
%         guidata(hObject,handles);
guidata(hObject,handles);
        drawnow;
end   
