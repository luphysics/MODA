%**************************************************************************
%***************************** Bispectrum GUI ******************************
%**************************************************************************
%---------------------------Credits----------------------------------------
% Wavelet Transform: Dmytro Iatsenko
% Wavelet bispectrum: Aleksandra Pidde
%----------------------------Documentation---------------------------------
% Coming Soon



function varargout = Bispectrum(varargin)
% BISPECTRUM MATLAB code for Bispectrum.fig
%      BISPECTRUM, by itself, creates a new BISPECTRUM or raises the existing
%      singleton*.
%
%      H = BISPECTRUM returns the handle to a new BISPECTRUM or the handle to
%      the existing singleton*.
%
%      BISPECTRUM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BISPECTRUM.M with the given input arguments.
%
%      BISPECTRUM('Property','Value',...) creates a new BISPECTRUM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bispectrum_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bispectrum_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help Bispectrum
% ------------------------------------------------------------------------
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


% Last Modified by GUIDE v2.5 12-Mar-2018 16:55:53
%*************************************************************************%
%                BEGIN initialization code - DO NOT EDIT                  %
%                ----------------------------------------                 %
%*************************************************************************%
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Bispectrum_OpeningFcn, ...
    'gui_OutputFcn',  @Bispectrum_OutputFcn, ...
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

function Bispectrum_OpeningFcn(hObject, eventdata, handles, varargin)
% Executes when GUI opens
handles=MODAsettings(hObject, handles);
handles.output = hObject;
set(handles.surr_plot,'Enable','off')
% set(handles.surr_num,'Enable','off')
% set(handles.alpha,'Enable','off')
guidata(hObject, handles);


function file_read_Callback(hObject, eventdata, handles)
% Load data into GUI
[handles,A]=MODAreadcheck(handles);
if A==1
    
    %-------Clearing axes and removing old data----------
    %Clearing all axes
    child_handles = allchild(handles.wt_pane);
    for i = 1:size(child_handles,1)
        if strcmp(get(child_handles(i),'type'),'axes')
            cla(child_handles(i),'reset');
            set(child_handles(i),'Visible','off');
        end
    end
    set(handles.display_type,'Enable','off');
    cla(handles.time_series_1,'reset')
    cla(handles.time_series_2,'reset')
    cla(handles.plot_pp,'reset')
    
    if isfield(handles, 'freqarr');handles = rmfield(handles, 'freqarr');    else    end
    if isfield(handles, 'sig');handles = rmfield(handles, 'sig');    else    end
    if isfield(handles, 'sig_cut');handles = rmfield(handles, 'sig_cut');    else    end
    if isfield(handles, 'biamp');handles = rmfield(handles, 'biamp');    else    end
    if isfield(handles, 'biamp');handles = rmfield(handles, 'biamp');    else    end
    if isfield(handles, 'biphase');handles = rmfield(handles, 'biphase');    else    end
    if isfield(handles, 'time_axis');handles = rmfield(handles, 'time_axis');    else    end
    if isfield(handles, 'pow_arr');handles = rmfield(handles, 'pow_arr');    else    end
    if isfield(handles, 'amp_arr');handles = rmfield(handles, 'amp_arr');    else    end
    if isfield(handles, 'pow_WT');handles = rmfield(handles, 'pow_WT');    else    end
    if isfield(handles, 'amp_WT');handles = rmfield(handles, 'amp_WT');    else    end
    if isfield(handles, 'WT');handles = rmfield(handles, 'WT');    else    end
    if isfield(handles, 'bispxxx');handles = rmfield(handles, 'bispxxx');    else    end
    if isfield(handles, 'bispxpp');handles = rmfield(handles, 'bispxpp');    else    end
    if isfield(handles, 'bisppxx');handles = rmfield(handles, 'bisppxx');    else    end
    if isfield(handles, 'bispppp');handles = rmfield(handles, 'bispppp');    else    end
    if isfield(handles, 'sampling_freq');handles = rmfield(handles, 'sampling_freq');    else    end
    if isfield(handles, 'peak_value');handles = rmfield(handles, 'peak_value');    else    end
    
    
    % Load data
    [handles]=MODAread(handles,1,"even");
    guidata(hObject,handles);
    
    if isfield(handles,'sig')
        
        % Plot time series
        plot(handles.time_series_1, handles.time_axis, handles.sig(1,:),'color',handles.linecol(1,:));%Plotting the time_series part afte calculation of appropriate limits
        xlim(handles.time_series_1,[0,size(handles.sig,2)./handles.sampling_freq]);
        ylabel(handles.time_series_1,'Sig 1','FontUnits','points','FontSize',10);
        
        if sum(abs(handles.sig(1,:)-handles.sig(2,:)))~=0
            plot(handles.time_series_2, handles.time_axis, handles.sig(2,:),'color',handles.linecol(1,:));%Plotting the time_series part afte calculation of appropriate limits
            xlim(handles.time_series_2,[0,size( handles.sig,2)./handles.sampling_freq]);
            
        else
        end
        
        xlabel(handles.time_series_2,'Time (s)','FontUnits','points','FontSize',10);
        ylabel(handles.time_series_2,'Sig 2','FontUnits','points','FontSize',10);
        
        set(handles.time_series_1,'FontUnits','points','FontSize',10);
        set(handles.time_series_1,'XTickLabels',[]);
        set(handles.time_series_2,'FontUnits','points','FontSize',10);
        linkaxes([handles.time_series_1 handles.time_series_2],'x');
        
        refresh_limits_Callback(handles.refresh_limits, eventdata, handles); % Updates data limits
        
        guidata(hObject,handles);
        preprocess_Callback(handles.preprocess, eventdata, handles);% Plots detrended curve
        set(handles.bisp_calc,'enable','on');
        guidata(hObject,handles);
        
        set(handles.status,'String','Data loaded. Continue with bispectrum calculation.');
    else
    end
else
    return;
end

function max_freq_Callback(hObject, eventdata, handles)
preprocess_Callback(hObject, eventdata, handles);

function min_freq_Callback(hObject, eventdata, handles)
preprocess_Callback(hObject, eventdata, handles);

function detrend_signal_popup_Callback(hObject, eventdata, handles)
% Detrends the signal and plots the chosen one
preprocess_Callback(handles.preprocess, eventdata, handles);

function handles=preprocess_Callback(hObject, eventdata, handles)
% Preprocesses and plots signal
ppstat=get(handles.preprocess,'Value');
if ppstat==2
    set(handles.plot_pp,'Visible','on')
    cla(handles.plot_pp,'reset');
    set(handles.detrend_signal_popup,'Enable','on')
    
    fmax = str2double(get(handles.max_freq,'String'));
    fmin = str2double(get(handles.min_freq,'String'));
    N=size(handles.sig_cut);
    L = N(2);
    handles.sig_pp=NaN(N);
    sig_select=get(handles.detrend_signal_popup,'Value');
    
    for j=1:2
        sig=handles.sig_cut(j,:);
        
        % Detrending
        X=(1:length(sig))'/handles.sampling_freq; XM=ones(length(X),4);
        
        for pn=1:3
            CX=X.^pn;
            XM(:,pn+1)=(CX-mean(CX))/std(CX);
        end
        sig=sig(:);
        new_sig=sig-XM*(pinv(XM)*sig);
        
        %Filtering
        fx=fft(new_sig,L); % Fourier transform of a signal
        
        Nq=ceil((L+1)/2);
        ff=[(0:Nq-1),-fliplr(1:L-Nq)]*handles.sampling_freq/L;
        ff=ff(:); % frequencies in Fourier transform
        
        fx(abs(ff)<=max([fmin,handles.sampling_freq/L]) | abs(ff)>=fmax)=0; % filter signal in a chosen frequency domain
        handles.sig_pp(j,:) = ifft(fx)';
        
    end
    % Plotting
    
    plot(handles.plot_pp,handles.time_axis_cut,handles.sig_cut(sig_select,:),'color',handles.linecol(1,:));
    hold(handles.plot_pp,'on');
    plot(handles.plot_pp,handles.time_axis_cut, handles.sig_pp(sig_select,:),'color',handles.linecol(2,:));
    legend(handles.plot_pp,{'Original','Pre-Processed'},'FontSize',8,'Location','Best','units','points');
    xlim(handles.plot_pp,[handles.time_axis_cut(1) handles.time_axis_cut(end)]);
    xlabel(handles.plot_pp,{'Time (s)'});
    
    drawnow;
    
else
    cla(handles.plot_pp,'reset');
    set(handles.plot_pp,'Visible','off')
    set(handles.detrend_signal_popup,'Enable','off')
end
guidata(hObject,handles);

function refresh_limits_Callback(hObject, eventdata, handles)
% Updates limits if user zooms in/out and presses refresh
x = get(handles.time_series_1,'xlim');
y = get(handles.time_series_1,'ylim');
xlim(handles.plot_pp,x)

t=x(2)-x(1);

xindex=x.*handles.sampling_freq;
x=strcat([num2str(x(1)),', ',num2str(x(2))]);
y=strcat([num2str(y(1)),', ',num2str(y(2))]);

xindex(2) = min(xindex(2),size(handles.sig,2));
xindex(1) = max(xindex(1),1);

handles.sig_cut=handles.sig(:,xindex(1):xindex(2));
handles.time_axis_cut = handles.time_axis(xindex(1):xindex(2));

set(handles.xlim,'String',x);
set(handles.ylim,'String',y);
set(handles.length,'String',t);

handles=preprocess_Callback(hObject, eventdata, handles);
handles.xl=[handles.time_axis_cut(1) handles.time_axis_cut(end)];
display_type_Callback(hObject, eventdata, handles)

guidata(hObject,handles);

function bisp_calc_Callback(hObject, eventdata, handles)
% Calculates the wavelet transform and the bispectrum
set(handles.bisp_calc,'Enable','off')

try
    fmax = str2double(get(handles.max_freq,'String'));
    fmin = str2double(get(handles.min_freq,'String'));
    fc =  str2double(get(handles.central_freq,'String'));
    fs = handles.sampling_freq;
    
    if isnan(fmax);fmax=fs/2;end
    
    % Preprocessing input
    x=get(handles.preprocess,'String'); ind=get(handles.preprocess,'Value'); ppselect=x{ind};
    
    if fmax>fs/2
        errordlg(['Maximum frequency cannot be higher than the Nyquist frequency. Please enter a value less than or equal to ',num2str(fs/2),' Hz.'],'Parameter Error');
        set(handles.bisp_calc,'Enable','on')
        return;
    end
    
    if fmin<=1/(length(handles.sig_cut)/handles.sampling_freq)
        errordlg('WT minimum frequency too low. To automatically calculate for minimum possible frequency leave "Min Freq" field blank.','Parameter Error');
        set(handles.bisp_calc,'Enable','on')
        return;
    end
    
    handles.status.String = "Testing";
    set(handles.status,'String','Calculating wavelet bispectrum');
    
    sig=handles.sig_cut;
    n = size(sig,1);
    handles.WT = cell(n, 1);
    handles.stop=0;
    
    handles.bispxxx=[];
    handles.bispppp=[];
    handles.bispxpp=[];
    handles.bisppxx=[];
    handles.surrxxx=[];
    handles.surrppp=[];
    handles.surrxpp=[];
    handles.surrpxx=[];
    
    sigcheck=sum(abs(handles.sig_cut(1,:)-handles.sig_cut(2,:)));
    handles.ns=str2num(get(handles.surr_num,'String'));
    
    if sigcheck~=0
        
        if strcmp(ppselect,'off')
            
            [handles.bispxxx] = bispecWavNew(handles.sig_cut(1,:),handles.sig_cut(1,:), fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',1,'wbar',1);
            handles.bispxxx=abs(handles.bispxxx);
            
            [handles.bispppp] = bispecWavNew(handles.sig_cut(2,:),handles.sig_cut(2,:), fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',2,'wbar',1);
            handles.bispppp=abs(handles.bispppp);
            
            [handles.bispxpp,handles.freqarr, handles.wavopt,WT1, WT2] = bispecWavNew(handles.sig_cut(1,:),handles.sig_cut(2,:), fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',3,'wbar',1);
            handles.bispxpp=abs(handles.bispxpp);
            
            [handles.bisppxx] = bispecWavNew(handles.sig_cut(2,:),handles.sig_cut(1,:), fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',4,'wbar',1);
            handles.bisppxx=abs(handles.bisppxx);
            
            if handles.ns>0
                
                set(handles.status,'String','Calculating wavelet bispectrum surrogates');
                
                handles.h = waitbar(0,'Calculating bispectrum surrogates...',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(handles.h,'canceling',0)
                guidata(hObject,handles);
                
                
                for j=1:handles.ns
                    waitbar(j / handles.ns,handles.h,sprintf(['Calculating surrogate (%d/%d)'],j,handles.ns));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    
                    surr1=wavsurrogate(handles.sig_cut(1,:),'IAAFT2',1);
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    
                    surr2=wavsurrogate(handles.sig_cut(2,:),'IAAFT2',1);
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    
                    handles.surrxxx(:,:,j)=abs(bispecWavNew(surr1,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',1));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    
                    handles.surrppp(:,:,j)=abs(bispecWavNew(surr2,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',2));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    
                    handles.surrxpp(:,:,j)=abs(bispecWavNew(surr1,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',3));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    
                    handles.surrpxx(:,:,j)=abs(bispecWavNew(surr2,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',4));
                    
                    
                end
                delete(handles.h);
            else
            end
        else
            [handles.bispxxx] = bispecWavNew(handles.sig_cut(1,:),handles.sig_cut(1,:),fs,'handles',handles,'hObject',hObject,'f0',fc,'fmin',fmin,'fmax',fmax,'num',1,'wbar',1);
            handles.bispxxx=abs(handles.bispxxx);
            
            [handles.bispppp] = bispecWavNew(handles.sig_cut(2,:),handles.sig_cut(2,:),fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',2,'wbar',1);
            handles.bispppp=abs(handles.bispppp);
            
            [handles.bispxpp,handles.freqarr,handles.wavopt,WT1, WT2] = bispecWavNew(handles.sig_cut(1,:),handles.sig_cut(2,:),fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',3,'wbar',1);
            handles.bispxpp=abs(handles.bispxpp);
            
            [handles.bisppxx] = bispecWavNew(handles.sig_cut(2,:),handles.sig_cut(1,:),fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',4,'wbar',1);
            handles.bisppxx=abs(handles.bisppxx);
            
            if handles.ns>0
                
                set(handles.status,'String','Calculating wavelet bispectrum surrogates');
                
                handles.h = waitbar(0,'Calculating bispectrum surrogates...',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(handles.h,'canceling',0)
                guidata(hObject,handles);
                for j=1:handles.ns
                    waitbar(j / handles.ns,handles.h,sprintf(['Calculating surrogate (%d/%d)'],j,handles.ns));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    surr1=wavsurrogate(handles.sig_cut(1,:),'IAAFT2',1);
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    surr2=wavsurrogate(handles.sig_cut(2,:),'IAAFT2',1);
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    handles.surrxxx(:,:,j)=abs(bispecWavNew(surr1,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',1));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    handles.surrppp(:,:,j)=abs(bispecWavNew(surr2,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',2));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    handles.surrxpp(:,:,j)=abs(bispecWavNew(surr1,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',3));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    handles.surrpxx(:,:,j)=abs(bispecWavNew(surr2,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',4));
                    
                    
                end
                delete(handles.h);
            else
            end
            
        end
    else
        linkaxes([handles.time_series_1, handles.time_series_2],'off')
        cla(handles.time_series_2,'reset')
        ylabel(handles.time_series_2,'Sig 2')
        ylabel(handles.time_series_1,'Sig 1')
        xlim(handles.time_series_2,[handles.time_axis_cut(1) handles.time_axis_cut(end)])
        linkaxes([handles.time_series_1, handles.time_series_2])
        box(handles.time_series_2,'on')
        if strcmp(ppselect,'off')
            [handles.bispxxx,handles.freqarr,handles.wavopt,WT1, WT2] = bispecWavNew(handles.sig_cut(1,:),handles.sig_cut(1,:), fs,'fmin',fmin,'fmax',fmax,'f0',fc,'preprocess','off','handles',handles,'hObject',hObject,'num',1,'wbar',1); handles.bispxxx=abs(handles.bispxxx);
            [handles.bispppp] = NaN(size(handles.bispxxx));
            [handles.bispxpp] = NaN(size(handles.bispxxx));
            [handles.bisppxx] = NaN(size(handles.bispxxx));
            
            if handles.ns>0
                set(handles.status,'String','Calculating wavelet bispectrum surrogates');
                
                handles.h = waitbar(0,'Calculating bispectrum surrogates...',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(handles.h,'canceling',0)
                guidata(hObject,handles);
                
                for j=1:handles.ns
                    waitbar(j / handles.ns,handles.h,sprintf(['Calculating surrogate (%d/%d)'],j,handles.ns));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    surr1=wavsurrogate(handles.sig_cut(1,:),'IAAFT2',1);
                    
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    handles.surrxxx(:,:,j)=abs(bispecWavNew(surr1,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',1));
                    handles.surrppp(:,:,j)=NaN(size(handles.surrxxx(:,:,j)));
                    handles.surrxpp(:,:,j)=NaN(size(handles.surrxxx(:,:,j)));
                    handles.surrpxx(:,:,j)=NaN(size(handles.surrxxx(:,:,j)));
                    
                end
                delete(handles.h);
            else
            end
        else
            [handles.bispxxx,handles.freqarr,handles.wavopt,WT1] = bispecWavNew(handles.sig_cut(1,:),handles.sig_cut(1,:),fs,'handles',handles,'hObject',hObject,'f0',fc,'fmin',fmin,'fmax',fmax,'num',1,'wbar',1); handles.bispxxx=abs(handles.bispxxx);
            WT2=NaN(size(WT1));
            [handles.bispppp] = NaN(size(handles.bispxxx));
            [handles.bispxpp] = NaN(size(handles.bispxxx));
            [handles.bisppxx] = NaN(size(handles.bispxxx));
            
            if handles.ns>0
                set(handles.status,'String','Calculating wavelet bispectrum surrogates');
                
                handles.h = waitbar(0,'Calculating bispectrum surrogates...',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(handles.h,'canceling',0)
                guidata(hObject,handles);
                
                for j=1:handles.ns
                    waitbar(j / handles.ns,handles.h,sprintf(['Calculating surrogate (%d/%d)'],j,handles.ns));
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    surr1=wavsurrogate(handles.sig_cut(1,:),'IAAFT2',1);
                    
                    if getappdata(handles.h,'canceling');guidata(hObject,handles);break;end
                    handles.surrxxx(:,:,j)=abs(bispecWavNew(surr1,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'handles',handles,'hObject',hObject,'num',1));
                    handles.surrppp(:,:,j)=NaN(size(handles.surrxxx(:,:,j)));
                    handles.surrxpp(:,:,j)=NaN(size(handles.surrxxx(:,:,j)));
                    handles.surrpxx(:,:,j)=NaN(size(handles.surrxxx(:,:,j)));
                end
                delete(handles.h);
            else
            end
        end
        
        
        
    end
    
    
    
    handles.amp_WT{1}=abs(WT1);
    handles.amp_arr{1}=nanmean(handles.amp_WT{1},2);
    handles.pow_WT{1}=abs(WT1).^2;
    handles.pow_arr{1}=nanmean(handles.pow_WT{1},2);
    
    handles.amp_WT{2}=abs(WT2);
    handles.amp_arr{2}=nanmean(handles.amp_WT{2},2);
    handles.pow_WT{2}=abs(WT2).^2;
    handles.pow_arr{2}=nanmean(handles.pow_WT{2},2);
    
    guidata(hObject,handles);
    drawnow;
    
    set(handles.display_type,'Enable','on')
    set(handles.save_mat,'Enable','on')
    set(handles.save_csv,'Enable','on')
    set(handles.save_session,'Enable','on')
    set(handles.bisp_calc,'Enable','on')
    set(handles.file_read,'Enable','off')
    set(handles.status,'String','Calculation complete');
    if handles.ns>0
        set(handles.surr_plot,'Enable','on')
    else
    end
    display_type_Callback(hObject, eventdata, handles)
    
catch e
    errordlg(e.message,'Error')
    set(handles.bisp_calc,'Enable','on')
    rethrow(e)
end

function display_type_Callback(hObject, eventdata, handles)

disp_select=get(handles.display_type,'Value');

set(handles.status,'String','Plotting data');

if length(handles.time_axis_cut)>=2000
    screensize = max(get(groot,'Screensize'));
    n = floor(size(handles.sig_cut,2)/screensize);
else
    n = 1;
end

%Clearing all axes
child_handles = allchild(handles.wt_pane);
for j = 1:size(child_handles,1)
    if strcmp(get(child_handles(j),'type'),'axes')
        cla(child_handles(j),'reset');
        set(child_handles(j),'Visible','off');
    end
end

% If wavelet transforms are selected
if (disp_select == 1 || disp_select == 2) && isfield(handles,'amp_WT')
    
    if length(handles.pow_WT{1})~=length(handles.time_axis_cut)
        errordlg('Time axis length has been changed, please recalculate bispectrum','Error')
        cla(handles.plot3d,'reset');
        cla(handles.plot_pow,'reset');
        set(handles.plot3d,'Visible','off');
        set(handles.plot_pow,'Visible','off');
        
    else
        set(handles.save_3dplot,'Enable','on')
        set(handles.save_both_plot,'Enable','on')
        set(handles.save_power_plot,'Enable','on')
        set(handles.save_bisp,'Enable','off')
        set(handles.save_biamp,'Enable','off')
        set(handles.save_biphase,'Enable','off')
        set(handles.save_bispect_biamp_biphase,'Enable','off')
        set(handles.save_all_plots,'Enable','off')
        set(handles.plot3d,'Visible','on');
        set(handles.plot_pow,'Visible','on');
        uistack(handles.plot3d,'top');
        uistack(handles.plot_pow,'top');
        
        
        if(handles.plot_type == 1)
            handles.peak_value = max(handles.pow_WT{disp_select,1}(:))+.1;
            pcolor(handles.plot3d, handles.time_axis_cut(1:n:end),handles.freqarr, handles.pow_WT{disp_select}(:,1:n:end));
            zlabel(handles.plot3d,'Power');
            xlabel(handles.plot_pow,'Average Power');
            plot(handles.plot_pow ,handles.pow_arr{disp_select}, handles.freqarr,'-k','LineWidth',3);
        else
            handles.peak_value = max(handles.amp_WT{disp_select}(:))+.1;
            pcolor(handles.plot3d, handles.time_axis_cut(1:n:end),handles.freqarr, handles.amp_WT{disp_select}(:,1:n:end));
            zlabel(handles.plot3d,'Amplitude');
            plot(handles.plot_pow ,handles.amp_arr{disp_select}, handles.freqarr,'-k','LineWidth',3);
            xlabel(handles.plot_pow,'Average Amplitude');
        end
        colormap(handles.plot3d,handles.cmap);
        shading(handles.plot3d,'interp');
        set(handles.plot3d,'yscale','log');
        set(handles.plot_pow,'yscale','log');
        set(handles.plot3d,'ylim',[min(handles.freqarr) max(handles.freqarr)]);
        set(handles.plot_pow,'ylim',[min(handles.freqarr) max(handles.freqarr)]);
        set(handles.plot3d,'xlim',handles.xl);
        xlabel(handles.plot3d,'Time (s)');
        ylabel(handles.plot3d,'Frequency (Hz)');
        ylabel(handles.plot_pow,'Frequency (Hz)');
        guidata(hObject,handles);
    end
    % If bispectra are selected
elseif disp_select > 2 && disp_select <7
    if length(handles.pow_WT{1})~=length(handles.time_axis_cut)
        errordlg('Time axis length has been changed, please recalculate bispectrum','Error')
        cla(handles.bisp,'reset');
        cla(handles.bisp_amp_axis,'reset');
        cla(handles.bisp_phase_axis,'reset');
        set(handles.bisp,'Visible','off');
        set(handles.bisp_amp_axis,'Visible','off');
        set(handles.bisp_phase_axis,'Visible','off');
        
        
    else
        set(handles.save_3dplot,'Enable','off')
        set(handles.save_both_plot,'Enable','off')
        set(handles.save_power_plot,'Enable','off')
        set(handles.save_bisp,'Enable','on')
        set(handles.save_biamp,'Enable','on')
        set(handles.save_biphase,'Enable','on')
        set(handles.save_bispect_biamp_biphase,'Enable','on')
        set(handles.save_all_plots,'Enable','off')
        cla(handles.bisp,'reset');
        cla(handles.bisp_amp_axis,'reset');
        cla(handles.bisp_phase_axis,'reset');
        
        set(handles.bisp,'Visible','on');
        set(handles.bisp_amp_axis,'Visible','on');
        set(handles.bisp_phase_axis,'Visible','on');
        uistack(handles.bisp,'top');
        uistack(handles.bisp_amp_axis,'top');
        uistack(handles.bisp_phase_axis,'top');
        surrselect=get(handles.surr_plot,'Value');
        
        if surrselect==1
            
            if disp_select == 3
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bispxxx)
                xlabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                title(handles.bisp,'Bispectrum 111');
            elseif disp_select == 4
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bispppp)
                xlabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                title(handles.bisp,'Bispectrum 222');
            elseif disp_select == 5
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bispxpp)
                xlabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                title(handles.bisp,'Bispectrum 122');
            elseif disp_select == 6
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bisppxx)
                xlabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                title(handles.bisp,'Bispectrum 211');
            end
        else
            
            handles.alph=str2num(get(handles.alpha,'String'));
            
            if floor((handles.ns+1)*handles.alph)==0
                K=1;
            else
                K=floor((handles.ns+1)*handles.alph);
            end
            
            if disp_select == 3
                S=sort(handles.surrxxx,3,'descend');
                
                handles.surrxxxT=S(:,:,K);
                handles.bispxxxS=handles.bispxxx-handles.surrxxxT;
                
                n=find(handles.bispxxxS<0);
                handles.bispxxxS(n)=NaN;
                
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bispxxxS)
                xlabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                title(handles.bisp,'Bispectrum 111');
            elseif disp_select == 4
                S=sort(handles.surrppp,3,'descend');
                handles.surrpppT=S(:,:,K);
                handles.bisppppS=handles.bispppp-handles.surrpppT;n=find(handles.bisppppS<0);handles.bisppppS(n)=NaN;
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bisppppS)
                xlabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                title(handles.bisp,'Bispectrum 222');
            elseif disp_select == 5
                S=sort(handles.surrxpp,3,'descend');
                handles.surrxppT=S(:,:,K);
                handles.bispxppS=handles.bispxpp-handles.surrxppT;n=find(handles.bispxppS<0);handles.bispxppS(n)=NaN;
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bispxppS)
                xlabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                title(handles.bisp,'Bispectrum 122');
            elseif disp_select == 6
                S=sort(handles.surrpxx,3,'descend');
                handles.surrpxxT=S(:,:,K);
                handles.bisppxxS=handles.bisppxx-handles.surrpxxT;n=find(handles.bisppxxS<0);handles.bisppxxS(n)=NaN;
                pcolor(handles.bisp, handles.freqarr, handles.freqarr, handles.bisppxxS)
                xlabel(handles.bisp,'Frequency - Sig 1 (Hz)');
                ylabel(handles.bisp,'Frequency - Sig 2 (Hz)');
                title(handles.bisp,'Bispectrum 211');
            end
        end
        
        set(handles.bisp, 'yscale','log','xscale','log');
        set(handles.bisp,'YTick',logspace(-5,5,11))
        set(handles.bisp,'XTick',logspace(-5,5,11))
        grid(handles.bisp,'on');
        
        idx_first = find(sum(~isnan(handles.bispxxx),1) > 0, 1 ,'first');
        idx_last = find(sum(~isnan(handles.bispxxx),1) > 0, 1 , 'last');
        xlim(handles.bisp,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
        ylim(handles.bisp,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
        
        colormap(handles.bisp,handles.cmap);
        shading(handles.bisp,'interp');
        set(handles.bisp,'fontunits','points','fontsize',10);
        
        
        
        guidata(hObject,handles);
    end
elseif disp_select == 7 && isfield(handles,'WT')
    set(handles.save_3dplot,'Enable','off')
    set(handles.save_both_plot,'Enable','off')
    set(handles.save_power_plot,'Enable','off')
    set(handles.save_bisp,'Enable','off')
    set(handles.save_biamp,'Enable','off')
    set(handles.save_biphase,'Enable','off')
    set(handles.save_bispect_biamp_biphase,'Enable','off')
    set(handles.save_all_plots,'Enable','on')
    set(handles.bispxxx_axis,'Visible','on');
    set(handles.bispppp_axis,'Visible','on');
    set(handles.bispxpp_axis,'Visible','on');
    set(handles.bisppxx_axis,'Visible','on');
    
    uistack(handles.bispxxx_axis,'top');
    uistack(handles.bispppp_axis,'top');
    uistack(handles.bispxpp_axis,'top');
    uistack(handles.bisppxx_axis,'top');
    surrselect=get(handles.surr_plot,'Value');
    if surrselect==1
        
        pcolor(handles.bispxxx_axis, handles.freqarr, handles.freqarr, handles.bispxxx)
        pcolor(handles.bispppp_axis, handles.freqarr, handles.freqarr, handles.bispppp)
        pcolor(handles.bisppxx_axis, handles.freqarr, handles.freqarr, handles.bisppxx)
        pcolor(handles.bispxpp_axis, handles.freqarr, handles.freqarr, handles.bispxpp)
    else
        handles.alph=str2num(get(handles.alpha,'String'));
        
        if floor((handles.ns+1)*handles.alph)==0
            K=1;
        else
            K=floor((handles.ns+1)*handles.alph);
        end
        S=sort(handles.surrxxx,3,'descend');
        handles.surrxxxT=S(:,:,K);
        handles.bispxxxS=handles.bispxxx-handles.surrxxxT;n=find(handles.bispxxxS<0);handles.bispxxxS(n)=NaN;
        pcolor(handles.bispxxx_axis, handles.freqarr, handles.freqarr, handles.bispxxxS)
        sigcheck=sum(abs(handles.sig_cut(1,:)-handles.sig_cut(2,:)));
        
        if sigcheck~=0
            
            S=sort(handles.surrppp,3,'descend');
            handles.surrpppT=S(:,:,K);
            handles.bisppppS=handles.bispppp-handles.surrpppT;n=find(handles.bisppppS<0);handles.bisppppS(n)=NaN;
            pcolor(handles.bispppp_axis, handles.freqarr, handles.freqarr, handles.bisppppS)
            
            S=sort(handles.surrxpp,3,'descend');
            handles.surrxppT=S(:,:,K);
            handles.bispxppS=handles.bispxpp-handles.surrxppT;n=find(handles.bispxppS<0);handles.bispxppS(n)=NaN;
            pcolor(handles.bispxpp_axis, handles.freqarr, handles.freqarr, handles.bispxppS)
            
            S=sort(handles.surrpxx,3,'descend');
            handles.surrpxxT=S(:,:,K);
            handles.bisppxxS=handles.bisppxx-handles.surrpxxT;n=find(handles.bisppxxS<0);handles.bisppxxS(n)=NaN;
            pcolor(handles.bisppxx_axis, handles.freqarr, handles.freqarr, handles.bisppxxS)
        else
        end
        
    end
    
    if length(handles.pow_WT{1})~=length(handles.time_axis_cut)
        errordlg('Time axis length has been changed, please recalculate bispectrum','Error')
        cla(handles.bispxxx_axis,'reset');
        cla(handles.bispppp_axis,'reset');
        cla(handles.bisppxx_axis,'reset');
        cla(handles.bispxpp_axis,'reset');
        set(handles.bispxxx_axis,'Visible','off');
        set(handles.bispppp_axis,'Visible','off');
        set(handles.bisppxx_axis,'Visible','off');
        set(handles.bispxpp_axis,'Visible','off');
    else
        
        if(handles.plot_type == 1)
            pcolor(handles.wt_1, handles.time_axis_cut(1:n:end) ,handles.freqarr, handles.pow_WT{1}(:,1:n:end));
            pcolor(handles.wt_2, handles.time_axis_cut(1:n:end) ,handles.freqarr, handles.pow_WT{2}(:,1:n:end));
        else
            pcolor(handles.wt_1, handles.time_axis_cut(1:n:end) ,handles.freqarr, handles.amp_WT{1}(:,1:n:end));
            pcolor(handles.wt_2, handles.time_axis_cut(1:n:end) ,handles.freqarr, handles.amp_WT{2}(:,1:n:end));
        end
        
        colormap(handles.wt_1,handles.cmap);
        colormap(handles.wt_2,handles.cmap);
        shading(handles.wt_1,'interp');
        shading(handles.wt_2,'interp');
        ylabel(handles.wt_1,'Frequency (Hz)');
        xlabel(handles.wt_2,'Time (s)');
        ylabel(handles.wt_2,'Frequency (Hz)');
        
        title(handles.wt_1,'Wavelet Transform - Signal 1');
        title(handles.wt_2,'Wavelet Transform - Signal 2');
        idx_first = find(sum(~isnan(handles.bispxxx),1) > 0, 1 ,'first');
        idx_last = find(sum(~isnan(handles.bispxxx),1) > 0, 1 , 'last');
        title(handles.bispxxx_axis,'b111');
        title(handles.bispppp_axis,'b222');
        title(handles.bisppxx_axis,'b211');
        title(handles.bispxpp_axis,'b122');
        ylabel(handles.bispxxx_axis,'Frequency(Hz)');
        ylabel(handles.bisppxx_axis,'Frequency(Hz)');
        xlabel(handles.bisppxx_axis,'Frequency(Hz)');
        xlabel(handles.bispppp_axis,'Frequency(Hz)');
        set(handles.wt_1,'yscale','log','ylim',[handles.freqarr(idx_first) handles.freqarr(idx_last)], 'xlim',handles.xl,...
            'zdir','reverse','xticklabel',[]);
        set(handles.wt_2,'yscale','log','ylim',[handles.freqarr(idx_first) handles.freqarr(idx_last)], 'xlim',handles.xl,...
            'zdir','reverse');
        
        child_handles = [handles.bispxxx_axis; handles.bispppp_axis; handles.bispxpp_axis; handles.bisppxx_axis];
        for i = 1:size(child_handles,1)
            if(strcmp(get(child_handles(i),'Type'),'axes'))
                colormap(child_handles(i),handles.cmap);
                shading(child_handles(i),'interp');
                set(child_handles(i),'yscale','log','xscale','log');
                
                xlim(child_handles(i),[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
                ylim(child_handles(i),[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
            end
        end
        
        guidata(hObject,handles);
    end
elseif ~isfield(handles,'bispxxx')
    return;
    
    
    
end
set(handles.status,'String','Plotting complete');

function mark_freq_Callback(hObject, eventdata, handles)
% Executes when pressing 'Select Point'
disp_select = get(handles.display_type,'Value');
if disp_select<=2 || disp_select>=7
    return;
end

[x,y] = ginput(1);

child_handles = allchild(handles.bisp);
for j = 1:size(child_handles,1)
    if(strcmp(get(child_handles(j),'Type'),'line'))
        xdat = get(child_handles(j),'XData');
        mark = get(child_handles(j),'Marker');
        if(length(xdat) == 1 && strcmp(mark,'*'))
            delete(child_handles(j))
        end
    end
end

hold(handles.bisp,'on');
plot(handles.bisp, x, y, 'k*')
set(handles.freq_1,'String',x);
set(handles.freq_2,'String',y);

function biph_calc_Callback(hObject, eventdata, handles)
% Executes on button press in biph_calc.
disp_select = get(handles.display_type,'Value');
if disp_select<=2 || disp_select>=7
    return;
end

if(disp_select>=3 && disp_select<=6)
    hold(handles.bisp_amp_axis,'on');
    hold(handles.bisp_phase_axis,'on');
    f1_str = get(handles.freq_1,'String');
    f2_str = get(handles.freq_2,'String');
    if (isempty(f1_str) || isempty(f2_str))
        return;
    end
    f1 = str2double(get(handles.freq_1,'String'));
    f2 = str2double(get(handles.freq_2,'String'));
    fmin=min(handles.freqarr);
    fmax=max(handles.freqarr);
    
    if f1<fmin || f1>fmax
        errordlg('Selected point is outside the allowable range','Parameter Error');
        return;
    end
    
    if f2<fmin || f2>fmax
        errordlg('Selected point is outside the allowable range','Parameter Error');
        return;
    end
    
    list = get(handles.frequency_select,'String');
    if isempty(list) %adds new cell to list containing f1 and f2
        list = cell(1,1);
        list{1,1} = sprintf(['%f, %f'],f1,f2);
    else
        list{size(list,1)+1,1} = sprintf(['%f, %f'],f1,f2);
    end
    sz = size(list,1);
    set(handles.frequency_select,'String',list,'Max',sz,'Value',sz);
    drawnow;
    
    
    %Marking the point
    child_handles = allchild(handles.bisp);
    for j = 1:size(child_handles,1)
        if(strcmp(get(child_handles(j),'Type'),'line'))
            xdat = get(child_handles(j),'XData');
            mark = get(child_handles(j),'Marker');
            if(length(xdat) == 1 && strcmp(mark,'*'))
                set(child_handles(j),'Marker','o','MarkerEdgeColor','k');
            end
        end
    end
    plot(handles.bisp,f1,f2,'ok');
    handles.freq_plot_list = {};
    for j = 1:size(list,1)
        temp = list{j,1};
        if temp(1,end) == 'n'
            handles.freq_plot_list{j,1} = temp(1:end-5);
        end
    end
    
    % Calculate biphase/biamplitude
    frequency_list = get(handles.frequency_select,'String');
    temp = cell2mat(frequency_list);
    F=str2num(temp);
    
    for j = 1:size(frequency_list,1)
        
        fl = F(j,:);
        fc =  str2double(get(handles.central_freq,'String'));
        % Preprocessing input
        x=get(handles.preprocess,'String'); ind=get(handles.preprocess,'Value'); ppselect=x{ind};
        if strcmp(ppselect,'off')
            
            if disp_select == 3
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_cut(1,:),handles.sig_cut(1,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            elseif disp_select == 4
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_cut(2,:),handles.sig_cut(2,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            elseif disp_select == 5
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_cut(1,:),handles.sig_cut(2,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            elseif disp_select == 6
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_cut(2,:),handles.sig_cut(1,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            end
        else
            if disp_select == 3
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_pp(1,:),handles.sig_pp(1,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            elseif disp_select == 4
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_pp(2,:),handles.sig_pp(2,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            elseif disp_select == 5
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_pp(1,:),handles.sig_pp(2,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            elseif disp_select == 6
                [handles.biamp{j}, handles.biphase{j}] = biphaseWavNew(handles.sig_pp(2,:),handles.sig_pp(1,:), handles.sampling_freq, [fl(2) fl(1)],handles.wavopt);
                xlabel(handles.bisp,'Frequency (Hz)');
                ylabel(handles.bisp,'Frequency (Hz)');
            end
        end
    end
    
    guidata(hObject, handles);
    frequency_select_Callback(hObject, eventdata, handles)
    
end

function plot_type_SelectionChangeFcn(hObject, eventdata, handles)
%deciding which plot
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'power'
        plot_type = 1;
    case 'amp'
        plot_type = 2;
end

handles.plot_type = plot_type;
guidata(hObject,handles);
display_selection = get(handles.display_type,'Value');
if display_selection<=2 || display_selection>=7
    display_type_Callback(handles.display_type, eventdata, handles)
    guidata(hObject,handles);
end

function frequency_select_Callback(hObject, eventdata, handles)
% Executes when selecting parameters in the frequency list
disp_select = get(handles.display_type,'Value');
if disp_select<=2 || disp_select>=7
    return;
elseif disp_select>=3 && disp_select<=6
    
    uistack(handles.bisp_amp_axis,'top');
    uistack(handles.bisp_phase_axis,'top');
    cla(handles.bisp_amp_axis);
    cla(handles.bisp_phase_axis);
    clear_axes_points(handles.bisp);
    hold(handles.bisp_phase_axis,'on');
    hold(handles.bisp_amp_axis,'on');
    hold(handles.bisp,'on');
    
    frequency_list = get(handles.frequency_select,'String');
    frequency_selected = get(handles.frequency_select,'Value');
    
    if isempty(frequency_list)
        return;
    end
    
    freq_selected_length = length(frequency_selected(1,:));
    handles.freq_plot_list = cell(freq_selected_length,1);
    handles.leg_bisp={};
    handles.f1list = {};
    handles.f2list = {};
    
    for j=1:freq_selected_length
        temp = cell2mat(frequency_list(frequency_selected(j),1));
        frequency_list{frequency_selected(j),1} = temp;
        handles.freq_plot_list{j,1} = temp;
        fl = csv_to_mvar(temp);
        
        colorindex = min(j,8); %so still plots lines even though no new linecolors available.
        plot(handles.bisp_amp_axis, handles.time_axis_cut, handles.biamp{frequency_selected(j)},'Linewidth',1,'color',handles.linecol(colorindex,:));
        grid(handles.bisp_amp_axis,'off');
        plot(handles.bisp_phase_axis, handles.time_axis_cut, handles.biphase{frequency_selected(j)},'Linewidth',1,'color',handles.linecol(colorindex,:));
        grid(handles.bisp_phase_axis,'off');
        ylabel(handles.bisp_amp_axis,'Biamplitude');
        ylabel(handles.bisp_phase_axis,'Biphase');
        xlabel(handles.bisp_phase_axis,'Time (s)');
        
        
        plot(handles.bisp, fl(1),fl(2), 'ok')
        
        shortened_list=get(handles.frequency_select,'String');
        shortened_list=cell2mat(shortened_list);
        
        handles.f1list{j}=str2double(shortened_list(frequency_selected(j),1:6));
        handles.f1list{j}=round(handles.f1list{j},3);
        handles.f1list{j}=mat2str(handles.f1list{j});
        
        handles.f2list{j}=str2double(shortened_list(frequency_selected(j),11:16));
        handles.f2list{j}=round(handles.f2list{j},3);
        handles.f2list{j}=mat2str(handles.f2list{j});
        handles.leg_bisp{j}=[handles.f1list{j},' - ',handles.f2list{j},' Hz'];
        
        
    end
    legend(handles.bisp_amp_axis,handles.leg_bisp,'FontSize',10);
    
end

set(handles.status,'String','Plotting complete');
guidata(hObject,handles);

function frequency_select_KeyPressFcn(hObject, eventdata, handles)
% Executes on key press in frequency_select
display_selection = get(handles.display_type,'Value');
if display_selection<=2 || display_selection>=7
    return;
end
switch eventdata.Key
    
    case 'delete'
        frequency_selected = get(handles.frequency_select,'Value');
        
        if min(frequency_selected)>1
            set(handles.frequency_select,'Value',min(frequency_selected)-1);
        else
            set(handles.frequency_select,'Value',1);
        end
        
        list = get(handles.frequency_select,'String');
        n=1:size(list,1);
        ne=n(1:end ~=frequency_selected);
        list=list(ne,:);
        set(handles.frequency_select,'String',list);
        drawnow;
        handles.biamp=handles.biamp(:,ne);
        handles.biphase=handles.biphase(:,ne);
        frequency_select_Callback(hObject, eventdata, handles)
end

function bisp_clear_Callback(hObject, eventdata, handles)
display_selection = get(handles.display_type,'Value');
if display_selection<=2 || display_selection>=7
    return;
end


cla(handles.bisp_amp_axis,'reset');
cla(handles.bisp_phase_axis,'reset');
set(handles.bisp_amp_axis,'fontunits','points','fontsize',10);
set(handles.bisp_phase_axis,'fontunits','points','fontsize',10);
clear_axes_points(handles.bisp);

%% Saving

function save_mat_Callback(hObject, eventdata, handles)
try
    [FileName,PathName] = uiputfile('.mat','Save data as');
    if FileName==0
        return;
    else
    end
    save_location = strcat(PathName,FileName)
    
    Bisp_data=savestrct(hObject, eventdata, handles);
    save(save_location,'Bisp_data');
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end

function save_csv_Callback(hObject, eventdata, handles)
try
    Bisp_data=savestrct(hObject, eventdata, handles);
    csvsavefolder(Bisp_data,handles)
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end

function csvsavefolder(D,handles)

curr=pwd;

[FileName,PathName] = uiputfile('.csv','Save as');
if FileName==0
    return;
else
end
cd(PathName)

foldername=[FileName(1:end-4)];
mkdir(foldername)



L=length(D.Time);



dstart=15;
data{1,1}='MODA v1.0 - Wavelet Bispectrum Analysis';
data{2,1}=date;
data{3,1}=[];
data{4,1}='PARAMETERS';
data{5,1}='Sampling frequency (Hz)';
data{5,2}=D.Sampling_frequency;
data{6,1}='Maximum frequency (Hz)';
data{6,2}=D.fmax;
data{7,1}='Minimum frequency (Hz)';
data{7,2}=D.fmin;
data{8,1}='Frequency resolution';
data{8,2}=D.fr;
data{9,1}='Preprocessing';
data{9,2}=D.Preprocessing;
data{10,1}='Cut Edges';
data{10,2}='on';
data{11,1}='Time start (s)';
data{11,2}=min(D.Time);
data{12,1}='Time end (s)';
data{12,2}=max(D.Time);



if isfield(D,'selected_plot')
    Np=size(D.selected_points,1);
    data{13,1}='Selected plot';
    data{13,2}=D.selected_plot;
    
    for j=1:Np
        data{13+j,1}=['Selected point ',num2str(j),' (Hz)'];
        data{13+j,2}=D.selected_points{j};
        dstart=dstart+1;
    end
    
    for k=1:Np
        data{dstart,k*2}=['Biamp - point ',num2str(k)];
        data{dstart,(k*2)+1}=['Biphase - point ',num2str(k)];
        for n=1:L
            data{n+dstart,k*2}=D.biamp{k}(n);
            data{n+dstart,(k*2)+1}=D.biphase{k}(n);
        end
    end
    data{dstart,1}='Time (s)';
    for l=1:L;
        data{l+dstart,1}=D.Time(l);
    end
    
    cell2csv([foldername,'\params_biphase_biamp.csv'],data,',');
else
    %         data{dstart,1}='Time (s)';
    %         for l=1:L;
    %             data{l+dstart,1}=D.Time(l);
    %         end
    %
    cell2csv([foldername,'\params.csv'],data,',');
end

sigcheck=sum(abs(handles.sig_cut(1,:)-handles.sig_cut(2,:)));

S=size(D.b111);
data2{1,1}='Freq';
for n=1:length(D.Frequency)
    data2{n+1,1}=D.Frequency(n);
    data2{1,n+1}=D.Frequency(n);
end
for j=1:S(1)
    for k=1:S(2)
        data2{j+1,k+1}=D.b111(j,k);
    end
end

cell2csv([foldername,'\Bispectrum_111.csv'],data2,',');
clear data2
if sigcheck ~=0
    data2{1,1}='Freq';
    for n=1:length(D.Frequency)
        data2{n+1,1}=D.Frequency(n);
        data2{1,n+1}=D.Frequency(n);
    end
    for j=1:S(1)
        for k=1:S(2)
            data2{j+1,k+1}=D.b222(j,k);
        end
    end
    
    cell2csv([foldername,'\Bispectrum_222.csv'],data2,',');
    clear data2
    
    data2{1,1}='Freq';
    for n=1:length(D.Frequency)
        data2{n+1,1}=D.Frequency(n);
        data2{1,n+1}=D.Frequency(n);
    end
    for j=1:S(1)
        for k=1:S(2)
            data2{j+1,k+1}=D.b122(j,k);
        end
    end
    
    cell2csv([foldername,'\Bispectrum_122.csv'],data2,',');
    clear data2
    
    data2{1,1}='Freq';
    for n=1:length(D.Frequency)
        data2{n+1,1}=D.Frequency(n);
        data2{1,n+1}=D.Frequency(n);
    end
    for j=1:S(1)
        for k=1:S(2)
            data2{j+1,k+1}=D.b211(j,k);
        end
    end
    cell2csv([foldername,'\Bispectrum_211.csv'],data2,',');
    
    clear data2
else
end

if isfield(D,'b111surr_threshold')
    data2{1,1}='Freq';
    for n=1:length(D.Frequency)
        data2{n+1,1}=D.Frequency(n);
        data2{1,n+1}=D.Frequency(n);
    end
    for j=1:S(1)
        for k=1:S(2)
            data2{j+1,k+1}=D.b111surr_threshold(j,k);
        end
    end
    cell2csv([foldername,'\Bispectrum_111_surr.csv'],data2,',');
    
    
    clear data2
    if sigcheck~=0
        data2{1,1}='Freq';
        for n=1:length(D.Frequency)
            data2{n+1,1}=D.Frequency(n);
            data2{1,n+1}=D.Frequency(n);
        end
        for j=1:S(1)
            for k=1:S(2)
                data2{j+1,k+1}=D.b222surr_threshold(j,k);
            end
        end
        cell2csv([foldername,'\Bispectrum_222_surr.csv'],data2,',');
        
        clear data2
        
        data2{1,1}='Freq';
        for n=1:length(D.Frequency)
            data2{n+1,1}=D.Frequency(n);
            data2{1,n+1}=D.Frequency(n);
        end
        for j=1:S(1)
            for k=1:S(2)
                data2{j+1,k+1}=D.b122surr_threshold(j,k);
            end
        end
        cell2csv([foldername,'\Bispectrum_122_surr.csv'],data2,',');
        
        clear data2
        
        data2{1,1}='Freq';
        for n=1:length(D.Frequency)
            data2{n+1,1}=D.Frequency(n);
            data2{1,n+1}=D.Frequency(n);
        end
        for j=1:S(1)
            for k=1:S(2)
                data2{j+1,k+1}=D.b211surr_threshold(j,k);
            end
        end
        cell2csv([foldername,'\Bispectrum_211_surr.csv'],data2,',');
        
    else
    end
else
end

cd(curr);

function Bisp_data=savestrct(hObject, eventdata, handles)
% Creates MATLAB structure for saving
% Preprocessing input
x=get(handles.preprocess,'String'); ind=get(handles.preprocess,'Value'); ppselect=x{ind};
sigcheck=sum(abs(handles.sig_cut(1,:)-handles.sig_cut(2,:)));
Bisp_data.b111=handles.bispxxx;
if sigcheck~=0
    Bisp_data.b222=handles.bispppp;
    Bisp_data.b122=handles.bispxpp;
    Bisp_data.b211=handles.bisppxx;
else
end

if isfield(handles,'alph')
    
    Bisp_data.b111surr=handles.surrxxx;
    if sigcheck~=0
        Bisp_data.b222surr=handles.surrppp;
        Bisp_data.b122surr=handles.surrxpp;
        Bisp_data.b211surr=handles.surrpxx;
    else
    end
    
    if floor((handles.ns+1)*handles.alph)==0
        K=1;
    else
        K=floor((handles.ns+1)*handles.alph);
    end
    
    S=sort(handles.surrxxx,3,'descend');
    handles.surrxxxT=S(:,:,K);
    if sigcheck ~=0
        S=sort(handles.surrppp,3,'descend');
        handles.surrpppT=S(:,:,K);
        S=sort(handles.surrxpp,3,'descend');
        handles.surrxppT=S(:,:,K);
        S=sort(handles.surrpxx,3,'descend');
        handles.surrpxxT=S(:,:,K);
    else
    end
    
    Bisp_data.b111surr_threshold=handles.surrxxxT;
    if sigcheck~=0
        Bisp_data.b222surr_threshold=handles.surrpppT;
        Bisp_data.b122surr_threshold=handles.surrxppT;
        Bisp_data.b211surr_threshold=handles.surrpxxT;
    else
    end
    Bisp_data.alpha=handles.alph;
    Bisp_data.surrnum=handles.ns;
    
else
end

frequency_list = get(handles.frequency_select,'String');

if ~isempty(frequency_list)
    Bisp_data.selected_points=frequency_list;
    display_selection = get(handles.display_type,'Value');
    
    if display_selection==3
        Bisp_data.selected_plot='111';
    elseif display_selection==4
        Bisp_data.selected_plot='222';
    elseif display_selection==5
        Bisp_data.selected_plot='122';
    elseif display_selection==6
        Bisp_data.selected_plot='211';
    end
    
    Bisp_data.biamp=handles.biamp;
    Bisp_data.biphase=handles.biphase;
else
end

if handles.plot_type==1
    avg_wt = cell2mat(handles.pow_arr);
    Bisp_data.WTPower=avg_wt';
else
    avg_wt = cell2mat(handles.amp_arr);
    Bisp_data.WTAmplitude=avg_wt';
end
Bisp_data.Frequency=handles.freqarr;
Bisp_data.Time=handles.time_axis_cut;
Bisp_data.Sampling_frequency=handles.sampling_freq;
Bisp_data.fmax=handles.freqarr(end);
Bisp_data.fmin=handles.freqarr(1);
Bisp_data.fr=str2double(get(handles.central_freq,'String'));
Bisp_data.Preprocessing=ppselect;


%% Plotting
function plot_TS_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.time_series_1, Fig);
ax2 = copyobj(handles.time_series_2, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.55,.85,.35]);
set(ax2,'Units', 'normalized', 'Position', [0.1,0.15,.85,.35]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);

function save_3dplot_Callback(hObject, eventdata, handles)
%Saves the 3d plot
Fig = figure;
ax = copyobj(handles.plot3d, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
colormap(ax,handles.cmap);
colorbar

function save_power_plot_Callback(hObject, eventdata, handles)
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
colorbar
ax2 = copyobj(handles.plot_pow, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.2,.55,.7]);
set(ax2,'Units', 'normalized', 'Position', [0.78,0.2,.18,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
ylabel(ax2,[])
set(Fig,'Units','normalized','Position', [0.2 0.2 0.6 0.5]);

colormap(Fig,handles.cmap);

function save_bisp_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.bisp, Fig);
set(ax,'Units', 'normalized', 'Position', [0.15,0.2,.7,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
set(Fig,'Units','normalized','Position', [0.2 0.2 0.4 0.5]);
colormap(ax,handles.cmap);
colorbar

function save_biamp_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.bisp_amp_axis, Fig);
set(ax,'Units', 'normalized', 'Position', [0.13,0.17,.8,.7], 'XTickMode','auto','YTickMode', 'auto','XTickLabelMode', 'auto', 'YTickLabelMode', 'auto');
xlabel(ax,'Time(s)');
set(Fig,'Units','normalized','Position', [0.3 0.3 0.5 0.45]);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
    set(lines(i),'LineWidth',2);
end
legend(ax,handles.leg_bisp,'FontSize',10);


function save_biphase_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.bisp_phase_axis, Fig);
set(ax,'Units', 'normalized', 'Position', [0.13,0.17,.8,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
set(Fig,'Units','normalized','Position', [0.3 0.3 0.5 0.45]);
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
    set(lines(i),'LineWidth',2);
end
legend(ax,handles.leg_bisp,'FontSize',10);


function save_bispect_biamp_biphase_Callback(hObject, eventdata, handles)
Fig = figure;

ax1 = copyobj(handles.bisp, Fig);
subplot(2,2,[1,3],ax1);
colorbar
colormap(Fig,handles.cmap);


ax2 = copyobj(handles.bisp_amp_axis, Fig);
subplot(2,2,2,ax2);
lines = findobj(gca,'Type','Line');
for i = 1:numel(lines)
    set(lines(i),'LineWidth',2);
end
set(ax2,'YTickMode', 'auto', 'YTickLabelMode', 'auto','XTickLabelMode','auto','FontUnits','points','FontSize',10);

ax3 = copyobj(handles.bisp_phase_axis, Fig);
subplot(2,2,4,ax3);
set(ax3,'YTickMode', 'auto', 'YTickLabelMode', 'auto');
lines = findobj(gca,'Type','Line');
for i = 1:numel(lines)
    set(lines(i),'LineWidth',2);
end
xlabel(ax3,'Time (s)','FontUnits','points','FontSize',10);
set(ax3,'YTickMode', 'auto', 'YTickLabelMode', 'auto','FontUnits','points','FontSize',10);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.9 0.5]);
legend(ax2,handles.leg_bisp,'FontSize',10);


function save_all_plots_Callback(hObject, eventdata, handles)
Fig = figure;

ax1 = copyobj(handles.bispxxx_axis, Fig);
subplot(2,3,1,ax1);
colormap(ax1,handles.cmap);

ax2 = copyobj(handles.bispxpp_axis, Fig);
subplot(2,3,2,ax2);
colormap(ax2,handles.cmap);

ax3 = copyobj(handles.bisppxx_axis, Fig);
subplot(2,3,4,ax3);
colormap(ax3,handles.cmap);

ax4 = copyobj(handles.bispppp_axis, Fig);
subplot(2,3,5,ax4);
colormap(ax4,handles.cmap);

ax5 = copyobj(handles.wt_1, Fig);
subplot(2,3,3,ax5);
colormap(ax5,handles.cmap);
set(ax5,'XTickLabelMode','auto','FontUnits','points','FontSize',10);

ax6 = copyobj(handles.wt_2, Fig);
subplot(2,3,6,ax6);
colormap(ax6,handles.cmap);

set(Fig,'Units','normalized','Position', [0.2 0.2 0.7 0.6]);

function varargout = Bispectrum_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function resetGUI_Callback(hObject, eventdata, handles)
% Executes when selecting 'New workspace'
Bispectrum;

function save_session_Callback(hObject, eventdata, handles)
% Executes when user selects 'Save current session'
MODAsave(handles)

function load_session_Callback(hObject, eventdata, handles)
% Executes when user selects 'Load previous session'
handles=MODAload;

function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Executes when user attempts to close GUI.
MODAclose(hObject,handles)

function surr_plot_Callback(hObject, eventdata, handles)
display_type_Callback(hObject, eventdata, handles)



function alpha_Callback(hObject, eventdata, handles)
display_type_Callback(hObject, eventdata, handles)
