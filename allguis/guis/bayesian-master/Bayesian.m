%Version 1.01 - Last edited by GL - 27/09/2017
%**************************************************************************
%***************************** Bayesian GUI ******************************
%**************************************************************************
%---------------------------Credits----------------------------------------
% Wavelet Transform: Dmytro Iatsenko
% Dynamical Bayesian Inference: Tomislav Stankovski
%----------------------------Documentation---------------------------------
%Coming Soon

function varargout = Bayesian(varargin)
% BAYESIAN MATLAB code for Bayesian.fig
%      BAYESIAN, by itself, creates a new BAYESIAN or raises the existing
%      singleton*.
%
%      H = BAYESIAN returns the handle to a new BAYESIAN or the handle to
%      the existing singleton*.
%
%      BAYESIAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BAYESIAN.M with the given input arguments.
%
%      BAYESIAN('Property','Value',...) creates a new BAYESIAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bayesian_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bayesian_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bayesian

% Last Modified by GUIDE v2.5 23-Jul-2018 15:43:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Bayesian_OpeningFcn, ...
    'gui_OutputFcn',  @Bayesian_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

function Bayesian_OpeningFcn(hObject,~,handles, varargin)
handles=MODAsettings(hObject, handles);
handles.c=0; % Increases with additional intervals added
handles.output = hObject;

% Disable plotting before data are loaded
set(handles.plot_TS,'Enable','off')
set(handles.plot_TSphi,'Enable','off')
set(handles.plot_phi1,'Enable','off')
set(handles.plot_phi2,'Enable','off')
set(handles.plot_CP,'Enable','off')
set(handles.plot_phiCS,'Enable','off')
set(handles.plot_CF,'Enable','off')
set(handles.save_csv,'Enable','off')
set(handles.save_mat,'Enable','off')
set(handles.cf_vid,'Enable','off')
set(handles.save_session,'Enable','off')
guidata(hObject, handles);

%----------------LOADING----------------------------------------------------

function file_read_Callback(hObject, eventdata, handles)
% Loads time series data
%% Load data
[handles]=MODAread(handles,1,"even");
if isfield(handles,'sampling_freq')
    refresh_limits_Callback(hObject, eventdata, handles);
else
end
guidata(hObject,handles); % Pass data overhead


function load_filt_Callback(hObject, eventdata, handles)
% Loads saved .mat file from filtering GUI, extracts the phases extracted
% from bands 1 and 2 and passes them to handles for use in the dynamical Bayesian
% inference
handles=MODAbayes_loadfilt(hObject,eventdata,handles,1);
if isfield(handles,'sampling_freq')
    refresh_limits_Callback(hObject, eventdata, handles);
else
end
guidata(hObject,handles); % Pass data overhead

function load_filt_2_Callback(hObject, eventdata, handles)
% Loads  2 saved .mat files from filtering GUI, extracts the phases extracted
% in eahc one and passes them to handles for use in the dynamical Bayesian
% inference
handles=MODAbayes_loadfilt(hObject,eventdata,handles,2);
if isfield(handles,'sampling_freq')
    refresh_limits_Callback(hObject, eventdata, handles);
else
end
guidata(hObject,handles); % Pass data overhead

function load_session_Callback(hObject, eventdata, handles)

handles=MODAload;

%---------------------------------------------------------------------------

function varargout = Bayesian_OutputFcn(~,~,handles)
varargout{1} = handles.output;

function handles=refresh_limits_Callback(hObject,eventdata,handles)
%% Updates limits and length in Graph Limits pane
x = get(handles.time_series_2,'xlim');
t = x(2) - x(1);
xstr = strcat([num2str(x(1)),', ',num2str(x(2))]);

xl = x.*handles.sampling_freq;

xl(2) = min(xl(2),size(handles.sig,2));
xl(1) = max(xl(1),1);

handles.sig_cut = handles.sig(:,xl(1):xl(2));

handles.time_axis_cut = handles.time_axis(:,xl(1):xl(2));

set(handles.xlim,'String',xstr);
set(handles.length,'String',t);
signal_list_Callback(hObject, eventdata, handles)

guidata(hObject,handles); % Pass data overhead

function signal_list_Callback(hObject, eventdata, handles)
% Selecting signal pair in signal list
sig_select = get(handles.signal_list, 'Value');
plot(handles.time_series_1, handles.time_axis_cut, handles.sig_cut(sig_select,:),'color',handles.linecol(1,:));
xl = csv_to_mvar(get(handles.xlim, 'String'));
xlim(handles.time_series_1, xl);
plot(handles.time_series_2, handles.time_axis_cut, handles.sig_cut(sig_select+size(handles.sig,1)/2,:),'color',handles.linecol(1,:));
xlim(handles.time_series_2, xl);
xlabel(handles.time_series_2, 'Time (s)');
ylabel(handles.time_series_1,'Sig 1');
ylabel(handles.time_series_2,'Sig 2');
set(handles.time_series_1,'xticklabel',[]);
%refresh_limits_Callback(hObject, eventdata, handles);
set(handles.status, 'String', 'Plot complete');

if isfield(handles,'p1') % Update all plots if the data exist
    %xyplot_Callback(hObject, eventdata, handles);
    interval_list_1_Callback(hObject, eventdata, handles)
else
    drawnow;
end


function add_interval_Callback(hObject,~,handles)
% Executes when 'Add interval' is pressed. Saves parameter set for
% calculations and displays info in the interval lists
handles.c=handles.c+1; % Increases index everytime interval is added
%% Prevents new data being loaded once parameters have been added (avoids bugs from overlapping data)
set(handles.file_read,'Enable','off')
set(handles.load_filt,'Enable','off')
set(handles.load_filt_2,'Enable','off')

%% Read frequency ranges and Bayesian parameters input by user
if isfield(handles,'pinput')
else
    f1r=csv_to_mvar(get(handles.freq_1,'String'));
    f2r=csv_to_mvar(get(handles.freq_2,'String'));
    
    if length(f1r)~=2
        errordlg('Incorrect format in Freq range 1. Example: to investigate 0.6 to 2 Hz, enter 0.6,2','Error');
        return;
    elseif length(f2r)~=2
        errordlg('Incorrect format in Freq range 2. Example: to investigate 0.6 to 2 Hz, enter 0.6,2','Error');
        return;
    else
    end
    
    handles.int1(handles.c,:) = f1r; % Frequency range 1
    handles.int2(handles.c,:) = f2r; % Frequency range 2
end
ws = str2double(get(handles.window_size,'String')); % Window size in seconds

% Default window size if no input from user (currently 10x longest period)
if isnan(ws)
    ws = 10/min([handles.int1(handles.c,:) handles.int2(handles.c,:)]);
end

handles.winds(handles.c,:)=ws;
handles.pr(handles.c,:)=str2double(get(handles.prop_const,'String')); % Propagation constant
handles.ovr(handles.c,:)=str2double(get(handles.overlap,'String')); % Overlap
handles.forder(handles.c,:)=str2double(get(handles.order,'String')); % Order of Fourier base function
handles.ns(handles.c,:)=str2double(get(handles.surrnum,'String')); % Number of surrogates
handles.confidence_level(handles.c,:)=str2double(get(handles.alphasig,'String'));

% handles.tm{1,handles.c}=NaN;
% handles.cc{1,handles.c}=NaN;
% handles.e{1,handles.c}=NaN;
% handles.cpl1{1,handles.c}=NaN;
% handles.cpl2{1,handles.c}=NaN;
% handles.cf1{1,handles.c}=NaN;
% handles.cf2{1,handles.c}=NaN;
% handles.mcf1{1,handles.c}=NaN;
% handles.mcf2{1,handles.c}=NaN;
% handles.surr_cpl1{1,handles.c}=NaN;
% handles.surr_cpl2{1,handles.c}=NaN;
% handles.p1{1,handles.c}=NaN;
% handles.p2{1,handles.c}=NaN;


%% Print selected frequency ranges and window size in GUI list

if isfield(handles,'pinput')
    fl = sprintf('%.3f,%.3f | %.2f | %.2f | %.2f | %d | %d | %d',min(str2num(handles.int1)),max(str2num(handles.int1)),ws,handles.ovr(handles.c,:),...
        handles.pr(handles.c,:),handles.forder(handles.c,:),handles.ns(handles.c,:),handles.confidence_level(handles.c,:));
    f2 = sprintf('%.3f,%.3f | %.2f | %.2f | %.2f | %d | %d | %d',min(str2num(handles.int2)),max(str2num(handles.int2)),ws,handles.ovr(handles.c,:),...
        handles.pr(handles.c,:),handles.forder(handles.c,:),handles.ns(handles.c,:),handles.confidence_level(handles.c,:));
else
    fl = sprintf('%.3f,%.3f | %.2f | %.2f | %.2f | %d | %d | %d',min(handles.int1(handles.c,:)),max(handles.int1(handles.c,:)),ws,handles.ovr(handles.c,:),...
        handles.pr(handles.c,:),handles.forder(handles.c,:),handles.ns(handles.c,:),handles.confidence_level(handles.c,:));
    f2 = sprintf('%.3f,%.3f | %.2f | %.2f | %.2f | %d | %d | %d',min(handles.int2(handles.c,:)),max(handles.int2(handles.c,:)),ws,handles.ovr(handles.c,:),...
        handles.pr(handles.c,:),handles.forder(handles.c,:),handles.ns(handles.c,:),handles.confidence_level(handles.c,:));
end

list = get(handles.interval_list_1,'String');
list{end+1,1} = fl;
set(handles.interval_list_1,'String',list);

list = get(handles.interval_list_2,'String');
list{end+1,1} = f2;
set(handles.interval_list_2,'String',list);


guidata(hObject, handles);

function calculate_Callback(hObject, eventdata, handles)
% Executes on 'Calculate'. Calculates dynamical Bayesian inference
x = get(handles.time_series_2,'xlim');
t = x(2) - x(1);
xstr = strcat([num2str(x(1)),', ',num2str(x(2))]);

xl = x.*handles.sampling_freq;

xl(2) = min(xl(2),size(handles.sig,2));
xl(1) = max(xl(1),1);

handles.sig_cut = handles.sig(:,xl(1):xl(2));

handles.time_axis_cut = handles.time_axis(:,xl(1):xl(2));
%% Prevents new data being loaded once calculation has been started (avoids bugs from overlapping data)
set(handles.file_read,'Enable','off')
set(handles.load_filt,'Enable','off')
set(handles.load_filt_2,'Enable','off')
set(handles.calculate,'Enable','off'); % User cannot press calculate during calculation

try
    set(handles.status,'String','Calculating...')
    drawnow;
    
    handles.h = waitbar(0,'Calculating DBI...',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(handles.h,'canceling',0)
    guidata(hObject,handles);
    N=size(handles.sig,1)/2;
    
    for k=1:handles.c % This is for the interval loop
        
        win=handles.winds(k);
        pr=handles.pr(k);
        ovr=handles.ovr(k);
        bn=handles.forder(k);
        ns=handles.ns(k);
        signif=handles.confidence_level(k);
        
        
        for j=1:size(handles.sig,1)/2
            if getappdata(handles.h,'canceling')
                break;
            end
            
            if isfield(handles,'pinput')
                phi1=handles.sig_cut(j,:);
                phi2=handles.sig_cut(j+N,:);
                handles.p1{j,k}=phi1;
                handles.p2{j,k}=phi2;
            else
                [handles.bands{j,k},~] = loop_butter(handles.sig_cut(j,:),handles.int1(k,:),handles.sampling_freq);
                phi1=angle(hilbert(handles.bands{j,k}));
                [handles.bands{j+N,k},~] = loop_butter(handles.sig_cut(j+N,:),handles.int2(k,:),handles.sampling_freq);
                phi2=angle(hilbert(handles.bands{j+N,k}));
                handles.p1{j,k}=phi1;
                handles.p2{j,k}=phi2;
                
            end
            
            
            %% Bayesian inference
            
            clear cpl1 cpl2 scpl1 scpl2
            [handles.tm{j,k},handles.cc{j,k},handles.e{j,k}]=bayes_main(phi1,phi2,win,1/handles.sampling_freq,ovr,pr,0,bn); % Main function
            for m=1:size(handles.cc{j,k},1)
                [cpl1(m),cpl2(m)]=dirc(handles.cc{j,k}(m,:),bn); % Direction of coupling
                [~,~,q21(:,:,m),q12(:,:,m)]=CFprint(handles.cc{j,k}(m,:),bn); % Coupling functions
            end
            
            handles.cpl1{j,k}=cpl1; % Coupling strength from 2 to 1
            handles.cpl2{j,k}=cpl2; % Coupling strength from 1 to 2
            handles.cf1{j,k}=q21; % Coupling functions for each time window
            handles.cf2{j,k}=q12;
            handles.mcf1{j,k}=squeeze(mean(q21,3)); % Mean coupling function
            handles.mcf2{j,k}=squeeze(mean(q12,3));
            
            
            
            %% Surrogates
            surr1=surrcalc(phi1,ns,'CPP',0,handles.sampling_freq);
            surr2=surrcalc(phi2,ns,'CPP',0,handles.sampling_freq);
            
            
            for n=1:ns
                [~,cc_surr{n}]=bayes_main(surr1(n,:),surr2(n,:),win,1/handles.sampling_freq,ovr,pr,1,bn);
                for idx=1:size(cc_surr{n},1)
                    [scpl1(n,idx),scpl2(n,idx)]=dirc(cc_surr{n}(idx,:),bn);
                end
            end
            
            %alph=str2num(get(handles.alphasig,'String'));
            alph=signif;
            alph=(1-(alph/100));
            if floor((ns+1)*alph)==0
                handles.surr_cpl1{j,k} = max(scpl1);
                handles.surr_cpl2{j,k} = max(scpl2);
            else
                K=floor((ns+1)*alph);
                s1=sort(scpl1,'descend');
                s2=sort(scpl2,'descend');
                handles.surr_cpl1{j,k} = s1(K,:);
                handles.surr_cpl2{j,k} = s2(K,:);
                
            end
            %         handles.surr_cpl1{j,k}=prctile(scpl1); % Define surrogate threshold here
            %         handles.surr_cpl2{j,k}=prctile(scpl2);
            
            waitbar((j+((size(handles.sig,1)/2)*(k-1)))/(size(handles.sig,1)/2*handles.c),handles.h,...
                sprintf(['Calculating DBI for signal pair ',num2str(j),', parameter set ',num2str(k)]));
            guidata(hObject,handles);
            
        end
        
    end
    delete(handles.h);
    guidata(hObject,handles);
    set(handles.status,'String','Calculation complete');
    set(handles.calculate,'Enable','on');
    set(handles.save_csv,'Enable','on')
    set(handles.save_mat,'Enable','on')
    set(handles.save_session,'Enable','on')
    interval_list_1_Callback(hObject, eventdata, handles)
catch e
    errordlg(e.message,'Error');
    set(handles.calculate,'Enable','on');
    delete(handles.h)
    rethrow(e)
    
    
end

function interval_list_1_Callback(hObject, eventdata, handles)

int_select = get(handles.interval_list_1,'Value');
set(handles.interval_list_2,'Value',int_select);
display_type_Callback(hObject, eventdata, handles)


function delete_set_Callback(hObject, eventdata, handles)
% --- Executes on button press in delete_set.
handles=MODAbayes_intdelete(hObject,eventdata,handles);
interval_list_1_Callback(hObject, eventdata, handles)


% Note: function is never called.
function interval_list_2_Callback(hObject, eventdata, handles)
int_select = get(handles.interval_list_2,'Value');
set(handles.interval_list_1,'Value',int_select);
display_type_Callback(hObject, eventdata, handles)


function display_type_Callback(hObject, eventdata, handles)
% Executes when display type drop down menu is changed
sig_select=get(handles.signal_list,'Value');
disp_select=get(handles.display_type,'Value');
int_select=get(handles.interval_list_1,'Value');
set(handles.interval_list_2,'Value',int_select);

% Remove dashed line markers
axes_child = allchild(handles.time_series_1);
axes_child2 = allchild(handles.time_series_2);
for j = 1:size(axes_child,1)
    if strcmpi(get(axes_child(j),'Type'),'Line')
        line_style = get(axes_child(j),'linestyle');
        if strcmp(line_style,'--')
            delete(axes_child(j));
            delete(axes_child2(j));
        end
    end
end

if disp_select==1
    % Enable plotting from menu of on-screen plots
    set(handles.plot_TSphi,'Enable','on')
    set(handles.plot_phi1,'Enable','on')
    set(handles.plot_phi2,'Enable','on')
    set(handles.plot_CP,'Enable','on')
    set(handles.plot_phiCS,'Enable','on')
    set(handles.plot_CF,'Enable','off')
    set(handles.curr_time,'visible','off');
    set(handles.cf_vid,'Enable','off')
    
    % Get the current x-limits of the time series plots to use later.
    xlim_backup = handles.time_series_1.XLim;
    
    % Clear axes
    child_handles = allchild(handles.plots_pane);
    for i = 1:length(child_handles)
        if strcmp(get(child_handles(i),'type'),'axes')
            cla(child_handles(i),'reset');
        end
    end
    
    % Set new x-limits of time series plots to their previous values,
    % to fix a bug where the signal is cut down to a tiny range after
    % using "delete parameter set" after a calculation has already been
    % attempted.
    handles.time_series_1.XLim = xlim_backup;
    
    set(handles.scaleon,'visible','off') % Remove match scale option
    set(handles.time_slider,'visible','off')
    set(handles.text23,'visible','off') % Time slider text
    set(handles.CF1,'visible','off');
    set(handles.CF2,'visible','off');
    set(handles.phi1_axes,'visible','on');
    set(handles.phi2_axes,'visible','on');
    set(handles.coupling_strength_axis,'visible','on');
    uistack(handles.phi1_axes,'top');
    uistack(handles.phi2_axes,'top');
    uistack(handles.coupling_strength_axis,'top');
    
    % Clear axes
    cla(handles.phi1_axes)
    cla(handles.phi2_axes)
    cla(handles.coupling_strength_axis)
    %linkaxes([handles.time_series_1 handles.time_series_2 handles.phi1_axes...
    %handles.phi2_axes handles.coupling_strength_axis],'x');
    
    
    
    % Removes dashed lines (may be leftover from plotting coupling functions)
    axes_child1=allchild(handles.time_series_1);
    axes_child2=allchild(handles.time_series_2);
    for j=1:size(axes_child1,1)
        if strcmpi(get(axes_child1(j),'Type'),'Line')
            line_style=get(axes_child1(j),'linestyle');
            if strcmp(line_style,'--')
                delete(axes_child1(j));
                delete(axes_child2(j));
            end
        end
    end
    if isfield(handles,'p1') && ~isempty(handles.p1)
        % Plot phases
        
        plot(handles.phi1_axes, handles.time_axis_cut, handles.p1{sig_select,int_select},'color',handles.linecol(1,:));
        plot(handles.phi2_axes, handles.time_axis_cut, handles.p2{sig_select,int_select},'color',handles.linecol(1,:));
        xlim(handles.phi1_axes, [handles.time_axis_cut(1) handles.time_axis_cut(end)]);
        xlim(handles.phi2_axes, [handles.time_axis_cut(1) handles.time_axis_cut(end)]);
        ylabel(handles.phi1_axes,'phi1');
        ylabel(handles.phi2_axes,'phi2');
        set(handles.phi1_axes,'xticklabel',[]);
        set(handles.phi2_axes,'xticklabel',[]);
        
        % Plot coupling strength
        hold(handles.coupling_strength_axis,'on');
        plot(handles.coupling_strength_axis,handles.tm{sig_select,int_select}+handles.time_axis_cut(1),handles.cpl1{sig_select,int_select},...
            'color',handles.linecol(1,:),'linewidth',handles.line2width);
        plot(handles.coupling_strength_axis,handles.tm{sig_select,int_select}+handles.time_axis_cut(1),handles.cpl2{sig_select,int_select},...
            'color',handles.linecol(2,:),'linewidth',handles.line2width);
        plot(handles.coupling_strength_axis,handles.tm{sig_select,int_select}+handles.time_axis_cut(1),handles.surr_cpl1{sig_select,int_select},...
            'color',handles.linecol(1,:),'linewidth',handles.line2width,'LineStyle','--');
        plot(handles.coupling_strength_axis,handles.tm{sig_select,int_select}+handles.time_axis_cut(1),handles.surr_cpl2{sig_select,int_select},...
            'color',handles.linecol(2,:),'linewidth',handles.line2width,'LineStyle','--');
        handles.leg1={'Data 2 -> 1','Data 1 -> 2','Surr 2 -> 1','Surr 1 -> 2'};
        legend(handles.coupling_strength_axis,handles.leg1,'orientation','horizontal');
        xlabel(handles.coupling_strength_axis,'Time (s)');
        ylabel(handles.coupling_strength_axis,'Coupling Strength');
        xlim(handles.coupling_strength_axis, [handles.time_axis_cut(1) handles.time_axis_cut(end)]);
        linkaxes([handles.time_series_1 handles.time_series_2 handles.phi1_axes handles.phi2_axes handles.coupling_strength_axis],'x'); % Ensures axis limits are identical for both plots
        
    else
    end
    
elseif disp_select==2
    % Disable plots and clear axes
    set(handles.plot_TSphi,'Enable','off')
    set(handles.plot_phi1,'Enable','off')
    set(handles.plot_phi2,'Enable','off')
    set(handles.plot_CP,'Enable','off')
    set(handles.plot_phiCS,'Enable','off')
    set(handles.plot_CF,'Enable','on')
    set(handles.cf_vid,'Enable','on')
    
    %     % Clear axes
    %     linkaxes([handles.time_series_1 handles.time_series_2 handles.phi1_axes...
    %         handles.phi2_axes handles.coupling_strength_axis],'off');
    
    cla(handles.phi1_axes,'reset');
    set(handles.phi1_axes,'visible','off');
    cla(handles.phi2_axes,'reset');
    set(handles.phi2_axes,'visible','off');
    cla(handles.coupling_strength_axis,'reset');
    set(handles.coupling_strength_axis,'visible','off');
    
    cla(handles.CF1);uistack(handles.CF1,'top');set(handles.CF1,'visible','on');
    cla(handles.CF2);uistack(handles.CF2,'top');set(handles.CF2,'visible','on');
    
    set(handles.scaleon,'visible','on')
    set(handles.time_slider,'visible','on')
    set(handles.text23,'visible','on') % Time slider text
    set(handles.time_slider,'Value',0)
    
    sig_select = get(handles.signal_list, 'Value');
    plot(handles.time_series_1, handles.time_axis_cut, handles.sig_cut(sig_select,:),'color',handles.linecol(1,:));
    xl = csv_to_mvar(get(handles.xlim, 'String'));
    xlim(handles.time_series_1, xl);
    plot(handles.time_series_2, handles.time_axis_cut, handles.sig_cut(sig_select+size(handles.sig,1)/2,:),'color',handles.linecol(1,:));
    xlim(handles.time_series_2, xl);
    xlabel(handles.time_series_2, 'Time (s)');
    ylabel(handles.time_series_1,'Sig 1');
    ylabel(handles.time_series_2,'Sig 2');
    set(handles.time_series_1,'xticklabel',[]);
    
    
    % Set up time slider
    
    minval=0;
    maxval=size(handles.tm{sig_select,int_select},2);
    stepsz=[1,5];
    set(handles.time_slider,'Max',maxval);
    set(handles.time_slider,'Min',minval);
    set(handles.time_slider,'SliderStep',stepsz/(maxval-minval));
    
    % Limits for surface plots
    t1=0:0.13:2*pi;t2=0:0.13:2*pi;
    cfmax1=max(max(handles.mcf1{sig_select,int_select}));
    cfmax2=max(max(handles.mcf2{sig_select,int_select}));
    
    if  cfmax1>cfmax2
        handles.cfmax=cfmax1;
    elseif cfmax2>cfmax1
        handles.cfmax=cfmax2;
    end
    
    cfmin1=min(min(handles.mcf1{sig_select,int_select}));
    cfmin2=min(min(handles.mcf2{sig_select,int_select}));
    
    if cfmin1<cfmin2
        handles.cfmin=cfmin1;
    elseif cfmin2<cfmin1
        handles.cfmin=cfmin2;
    end
    
    % Plot coupling functions
    set(handles.curr_time,'visible','off'); % Turn off current time for mean plots
    surf(handles.CF1,t1,t2,handles.mcf1{sig_select,int_select},'FaceColor','interp');
    xlabel(handles.CF1,'\phi_1');ylabel(handles.CF1,'\phi_2');zlabel(handles.CF1,'q_1(\phi_1,\phi_2)');
    title(handles.CF1,'Time-averaged CF 2 --> 1');
    view(handles.CF1,[-40 50]);
    colormap(handles.CF1,handles.cmap)
    
    surf(handles.CF2,t1,t2,handles.mcf2{sig_select,int_select},'FaceColor','interp');
    xlabel(handles.CF2,'\phi_1');ylabel(handles.CF2,'\phi_2');zlabel(handles.CF2,'q_2(\phi_1,\phi_2)');
    title(handles.CF2,'Time-averaged CF 1 --> 2');
    view(handles.CF2,[-40 50]);
    colormap(handles.CF2,handles.cmap)
    
    xlim(handles.CF1,[0 2*pi]);
    ylim(handles.CF1,[0 2*pi]);
    xlim(handles.CF2,[0 2*pi]);
    ylim(handles.CF2,[0 2*pi]);
    
    if get(handles.scaleon,'Value')==1 % Z axis scaling
        set(handles.CF1,'Zlim',[handles.cfmin handles.cfmax])
        set(handles.CF2,'Zlim',[handles.cfmin handles.cfmax])
    else
    end
    
    
else
end
guidata(hObject,handles);


function scaleon_Callback(hObject, eventdata, handles)
% --- Executes on button press in scaleon.

time_slider_Callback(hObject, eventdata, handles)

function time_slider_Callback(hObject, eventdata, handles)

disp_select=get(handles.display_type,'Value');

if disp_select~=2
    return;
else
end

sig_select=get(handles.signal_list,'Value');
int_select=get(handles.interval_list_1,'Value');
win=handles.winds(int_select);


slider_val=round(get(handles.time_slider,'Value'));
if slider_val==0
    % Remove dashed line markers
    axes_child = allchild(handles.time_series_1);
    axes_child2 = allchild(handles.time_series_2);
    for j = 1:size(axes_child,1)
        if strcmpi(get(axes_child(j),'Type'),'Line')
            line_style = get(axes_child(j),'linestyle');
            if strcmp(line_style,'--')
                delete(axes_child(j));
                delete(axes_child2(j));
            end
        end
    end
    interval_list_1_Callback(hObject, eventdata, handles)
    return;
else
end

t1=0:0.13:2*pi;t2=0:0.13:2*pi;
q1(1:length(t1),1:length(t1))=0;q2=q1;
u = handles.cc{sig_select,int_select}(slider_val,:);
K=length(u)/2;

% Remove dashed line markers
axes_child = allchild(handles.time_series_1);
axes_child2 = allchild(handles.time_series_2);
for j = 1:size(axes_child,1)
    if strcmpi(get(axes_child(j),'Type'),'Line')
        line_style = get(axes_child(j),'linestyle');
        if strcmp(line_style,'--')
            delete(axes_child(j));
            delete(axes_child2(j));
        end
    end
end

% Update dashed line markers
yl = get(handles.time_series_1,'ylim');
y2 = get(handles.time_series_2,'ylim');
handles.tp1=(handles.tm{sig_select,int_select}(slider_val)-win/2)+handles.time_axis_cut(1);
xl=[handles.tp1 handles.tp1];
hold(handles.time_series_1,'on');
plot(handles.time_series_1,xl,yl,'--','color',[0.8500    0.3250    0.0980],'linewidth',1)
hold(handles.time_series_1,'off');
hold(handles.time_series_2,'on');
plot(handles.time_series_2,xl,y2,'--','color',[0.8500    0.3250    0.0980],'linewidth',1)
hold(handles.time_series_2,'off');

yl = get(handles.time_series_1,'ylim');
y2 = get(handles.time_series_2,'ylim');
handles.tp2=(handles.tm{sig_select,int_select}(slider_val)+win/2)+handles.time_axis_cut(1);
xl=[handles.tp2 handles.tp2];
hold(handles.time_series_1,'on');
plot(handles.time_series_1,xl,yl,'--','color',[0.8500    0.3250    0.0980],'linewidth',1)
hold(handles.time_series_1,'off');
hold(handles.time_series_2,'on');
plot(handles.time_series_2,xl,y2,'--','color',[0.8500    0.3250    0.0980],'linewidth',1)
hold(handles.time_series_2,'off');

set(handles.scaleon,'visible','on');
cfmax1=max(max(handles.cf1{sig_select,int_select}(:,:,slider_val)));
cfmax2=max(max(handles.cf2{sig_select,int_select}(:,:,slider_val)));

if  cfmax1>cfmax2
    handles.cfmax=cfmax1;
elseif cfmax2>cfmax1
    handles.cfmax=cfmax2;
end

cfmin1=min(min(handles.cf1{sig_select,int_select}(:,:,slider_val)));
cfmin2=min(min(handles.cf2{sig_select,int_select}(:,:,slider_val)));

if cfmin1<cfmin2
    handles.cfmin=cfmin1;
elseif cfmin2<cfmin1
    handles.cfmin=cfmin2;
end

surf(handles.CF1,t1,t2,handles.cf1{sig_select,int_select}(:,:,slider_val),'FaceColor','interp');
surf(handles.CF2,t1,t2,handles.cf2{sig_select,int_select}(:,:,slider_val),'FaceColor','interp');
xlabel(handles.CF1,'\phi_1');ylabel(handles.CF1,'\phi_2');zlabel(handles.CF1,'q_1(\phi_1,\phi_2)');
xlabel(handles.CF2,'\phi_1');ylabel(handles.CF2,'\phi_2');zlabel(handles.CF2,'q_2(\phi_1,\phi_2)');
title(handles.CF1,'CF 2 --> 1')
title(handles.CF2,'CF 1 --> 2')
xlim(handles.CF1,[0 2*pi]);
ylim(handles.CF1,[0 2*pi]);
xlim(handles.CF2,[0 2*pi]);
ylim(handles.CF2,[0 2*pi]);

if get(handles.scaleon,'Value')==1
    set(handles.CF1,'Zlim',[handles.cfmin handles.cfmax])
    set(handles.CF2,'Zlim',[handles.cfmin handles.cfmax])
else
end
set(handles.curr_time,'visible','on');
set(handles.curr_time,'string',[num2str(handles.tp1),' - ',num2str(handles.tp2),' s']);
view(handles.CF1,[-40 50]);
view(handles.CF2,[-40 50]);
guidata(hObject,handles);


%% Saving callbacks

function save_csv_Callback(hObject, eventdata, handles)
try
    Bayes_data.phase1=handles.p1;
    Bayes_data.phase2=handles.p2;
    Bayes_data.sampling_freq=handles.sampling_freq;
    Bayes_data.time=handles.time_axis_cut;
    Bayes_data.interval1=handles.int1;
    Bayes_data.interval2=handles.int2;
    Bayes_data.surrnum=handles.ns;
    Bayes_data.confidence_level=handles.confidence_level;
    Bayes_data.Bayes_win=handles.winds;
    Bayes_data.overlap=handles.ovr;
    Bayes_data.propagation_const=handles.pr;
    Bayes_data.Fourier_base=handles.forder;
    Bayes_data.Bayestime=handles.tm;
    % Bayes_data.Bayes_inferred=handles.cc;
    % Bayes_data.Bayes_noise=handles.e;
    Bayes_data.coupling_strength_2to1=handles.cpl1;
    Bayes_data.coupling_strength_1to2=handles.cpl2;
    Bayes_data.coupling_function_2to1=handles.cf1;
    Bayes_data.coupling_function_1to2=handles.cf2;
    Bayes_data.mean_cfunc_2to1=handles.mcf1;
    Bayes_data.mean_cfunc_1to2=handles.mcf2;
    Bayes_data.surrogate_coupling_strength_2to1=handles.surr_cpl1;
    Bayes_data.surrogate_coupling_strength_1to2=handles.surr_cpl2;
    
    csvsavefolder(Bayes_data)
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end


% --------------------------------------------------------------------
function save_mat_Callback(hObject, eventdata, handles)
try
    [FileName,PathName] = uiputfile('.mat','Save as');
    if FileName==0
        return;
    else
    end
    save_location = strcat(PathName,FileName);
    
    Bayes_data.phase1=handles.p1;
    Bayes_data.phase2=handles.p2;
    Bayes_data.sampling_freq=handles.sampling_freq;
    Bayes_data.time=handles.time_axis_cut;
    Bayes_data.interval1=handles.int1;
    Bayes_data.interval2=handles.int2;
    Bayes_data.surrnum=handles.ns;
    Bayes_data.confidence_level=handles.confidence_level;
    Bayes_data.Bayes_win=handles.winds;
    Bayes_data.overlap=handles.ovr;
    Bayes_data.propagation_const=handles.pr;
    Bayes_data.Fourier_base=handles.forder;
    Bayes_data.Bayestime=handles.tm;
    % Bayes_data.Bayes_inferred=handles.cc;
    % Bayes_data.Bayes_noise=handles.e;
    Bayes_data.coupling_strength_2to1=handles.cpl1;
    Bayes_data.coupling_strength_1to2=handles.cpl2;
    Bayes_data.coupling_function_2to1=handles.cf1;
    Bayes_data.coupling_function_1to2=handles.cf2;
    Bayes_data.mean_cfunc_2to1=handles.mcf1;
    Bayes_data.mean_cfunc_1to2=handles.mcf2;
    Bayes_data.surrogate_coupling_strength_2to1=handles.surr_cpl1;
    Bayes_data.surrogate_coupling_strength_1to2=handles.surr_cpl2;
    
    save(save_location,'Bayes_data');
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end

function cf_vid_Callback(hObject, eventdata, handles)
% Saves video of coupling functions currently selected
int_select = get(handles.interval_list_1,'Value');
sig_select = get(handles.signal_list,'Value');
c1=handles.cf1{sig_select,int_select};
c2=handles.cf2{sig_select,int_select};
win=handles.winds(int_select);
t=handles.tm{sig_select,int_select};



[FileName,PathName] = uiputfile('.avi','Save video as');
save_location = strcat(PathName,FileName)

if isempty(save_location)
    return;
else
end

v = VideoWriter(save_location);
v.FrameRate=2;
open(v)

figure('position',[100 100 1200 500])

for j=1:size(c1,3)
    
    cfmax1=max(max(c1(:,:,j)));
    cfmax2=max(max(c2(:,:,j)));
    
    if  cfmax1>cfmax2
        handles.cfmax=cfmax1;
    elseif cfmax2>cfmax1
        handles.cfmax=cfmax2;
    end
    
    cfmin1=min(min(c1(:,:,j)));
    cfmin2=min(min(c2(:,:,j)));
    
    if cfmin1<cfmin2
        handles.cfmin=cfmin1;
    elseif cfmin2<cfmin1
        handles.cfmin=cfmin2;
    end
    
    
    subplot(1,2,1)
    surf(squeeze(c1(:,:,j)))
    title('Oscillator 2 to oscillator 1')
    view([-40 50])
    set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
    xlabel('\phi_1');ylabel('\phi_2');zlabel('q_1(\phi_1,\phi_2)');axis tight
    grid on
    if get(handles.scaleon,'Value')==1
        set(gca,'Zlim',[handles.cfmin handles.cfmax])
    else
    end
    subplot(1,2,2)
    surf(squeeze(c2(:,:,j)))
    title('Oscillator 1 to oscillator 2')
    view([-40 50])
    set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
    xlabel('\phi_1');ylabel('\phi_2');zlabel('q_2(\phi_1,\phi_2)');axis tight
    grid on
    if get(handles.scaleon,'Value')==1
        set(gca,'Zlim',[handles.cfmin handles.cfmax])
    else
    end
    colormap(handles.cmap)
    set(gcf,'renderer','zbuffer')
    
    ax2=axes('Position',[.45 .1 .15 .15],'visible','off');
    %currt=[num2str(t{sig_select,int_select}-floor(win/2)),' - ',num2str(t{sig_select,int_select}+floor(win/2)),' s'];
    currt=[num2str(t(j)-floor(win/2)),' - ',num2str(t(j)+floor(win/2)),' s'];
    text(0.2, 0.2,currt,'fontsize',12,'fontweight','bold')
    set(gca,'nextplot','replacechildren');
    currFrame = getframe(gcf);
    writeVideo(v,currFrame)
    
end
close(v)
close(gcf)

function csvsavefolder(Bayes_data)
% Saving function
curr=pwd;

[FileName,PathName] = uiputfile('.csv','Save as');
if FileName==0
    return;
else
end
cd(PathName)

foldername=[FileName(1:end-4)];
mkdir(foldername)



Nb=size(Bayes_data.interval1,1);

for nn=1:Nb
    comp=nn;
    clear data
    Ns=size(Bayes_data.phase1,1);
    L=length(Bayes_data.Bayestime{1,comp});
    
    
    
    data{1,1}='MODA v1.0 - Dynamical Bayesian Inference';
    data{2,1}=date;
    data{3,1}=[];
    data{4,1}='PARAMETERS';
    data{5,1}='Sampling frequency (Hz)';
    data{5,2}=Bayes_data.sampling_freq;
    data{6,1}='Time start (s)';
    data{6,2}=min(Bayes_data.time);
    data{7,1}='Time end (s)';
    data{7,2}=max(Bayes_data.time);
    
    for n=1
        
        data{8,1}='Overlap';
        data{8,2}=Bayes_data.overlap(comp);
        data{9,1}='Propagation constant';
        data{9,2}=Bayes_data.propagation_const(comp);
        data{10,1}='Fourier basis';
        data{10,2}=Bayes_data.Fourier_base(comp);
        data{11,1}='Number of surrogates';
        data{11,2}=Bayes_data.surrnum(comp);
        data{12,1}='Window size (s)';
        data{12,n+1}=Bayes_data.Bayes_win(comp);
        data{13,1}=['Interval 1'];
        data{13,n+1}=Bayes_data.interval1(comp,:);
        data{14,1}=['Interval 2'];
        data{14,n+1}=Bayes_data.interval2(comp,:);
        data{15,1}='Confidence level';
        data{15,n+1}=Bayes_data.confidence_level(comp,:);
        c=n+15;
    end
    
    dstart=c+1;
    
    data{dstart,1}='Time (s)';
    for l=1:L;
        data{l+dstart,1}=Bayes_data.Bayestime{1,comp}(l);
    end
    
    for h=1:Ns
        data{dstart,h*2}=['Subject ',num2str(h),' - 2 to 1'];
        data{dstart,h*2+1}=['Subject ',num2str(h),' - 1 to 2'];
        data{dstart,h*2+Ns*2}=['Surrogate ',num2str(h),' - 2 to 1'];
        data{dstart,h*2+Ns*2+1}=['Surrogate ',num2str(h),' - 1 to 2'];
        for l=1:L;
            data{l+dstart,h*2}=Bayes_data.coupling_strength_2to1{h,comp}(l);
            data{l+dstart,h*2+1}=Bayes_data.coupling_strength_1to2{h,comp}(l);
            data{l+dstart,h*2+Ns*2}=Bayes_data.surrogate_coupling_strength_2to1{h,comp}(l);
            data{l+dstart,h*2+Ns*2+1}=Bayes_data.surrogate_coupling_strength_1to2{h,comp}(l);
        end
    end
    
    cell2csv([foldername,'\Coupling_strength',num2str(nn),'.csv'],data,',');
    
end

t1=0:0.13:2*pi;t2=0:0.13:2*pi;

for b=1:Ns
    for bb=1:Nb
        meanCF1=Bayes_data.mean_cfunc_2to1{b,bb};
        meanCF2=Bayes_data.mean_cfunc_1to2{b,bb};
        %         S=size(Bayes_data.mean_cfunc_2to1{1,1},1);
        %
        %         M(2:S+1,2:S+1)=meanCF1;
        %         M(2:S+1,2+S+1:S*2+2)=meanCF2;
        %         M(2:S+1,1)=t1;
        %         M(1,2:S+1)=t2;
        %         M(1,2+S+1:S*2+2)=t2;
        %         csvwrite([foldername,'\Mean_CF ',num2str(bb),' - Subject ',num2str(b),'.csv'],M)
        %
        %     end
        % end
        
        S=size(meanCF1);
        for n=1:length(t1)
            M{n+1,1}=t1(n);
            M{1,n+1}=t2(n);
            M{1,n+2+length(t1)}=t2(n);
        end
        
        for j=1:S(1)
            for k=1:S(2)
                M{j+1,k+1}=meanCF1(j,k);
                M{j+1,k+2+length(t1)}=meanCF2(j,k);
            end
        end
        cell2csv([foldername,'\Mean_CF ',num2str(bb),' - Subject ',num2str(b),'.csv'],M,',')
    end
end

cd(curr);

function save_session_Callback(hObject, eventdata, handles)

MODAsave(handles)

%% Plotting callbacks

function plot_TS_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.time_series_1, Fig);
ax2 = copyobj(handles.time_series_2, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.55,.85,.35]);
set(ax2,'Units', 'normalized', 'Position', [0.1,0.15,.85,.35]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);

% --------------------------------------------------------------------
function plot_TSphi_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.time_series_1, Fig);
ax2 = copyobj(handles.phi1_axes, Fig);
ax3 = copyobj(handles.time_series_2, Fig);
ax4 = copyobj(handles.phi2_axes, Fig);
xlabel('Time (s)')
set(ax4, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
set(ax3, 'XTickLabel', [])
xlabel(ax3,'')
set(ax1,'Units', 'normalized', 'Position', [0.1,0.79,.85,.17]);
set(ax2,'Units', 'normalized', 'Position', [0.1,0.56,.85,.17]);
set(ax3,'Units', 'normalized', 'Position', [0.1,0.33,.85,.17]);
set(ax4,'Units', 'normalized', 'Position', [0.1,0.1,.85,.17]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);

% --------------------------------------------------------------------
function plot_phi1_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.phi1_axes, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(ax1, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
xlabel('Time (s)')
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.25]);


% --------------------------------------------------------------------
function plot_phi2_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.phi2_axes, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(ax1, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
xlabel('Time (s)')
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.25]);


% --------------------------------------------------------------------
function plot_CP_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.coupling_strength_axis, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(ax1, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
xlabel('Time (s)')
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.25]);
legend(handles.leg1)

% --------------------------------------------------------------------
function plot_phiCS_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.phi1_axes, Fig);
ax2 = copyobj(handles.phi2_axes, Fig);
ax3 = copyobj(handles.coupling_strength_axis, Fig);
xlabel('Time (s)')
set(ax1,'Units', 'normalized', 'Position', [0.1,0.7,.85,.27]);
set(ax2,'Units', 'normalized', 'Position', [0.1,0.4,.85,.27]);
set(ax3,'Units', 'normalized', 'Position', [0.1,0.1,.85,.27]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.4]);
legend(handles.leg1,'orientation','horizontal')


% --------------------------------------------------------------------
function plot_CF_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.CF1, Fig);
ax2 = copyobj(handles.CF2, Fig);
ax3 = copyobj(handles.curr_time, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.17,.35,.75]);
set(ax2,'Units', 'normalized', 'Position', [0.6,0.17,.35,.75]);
set(ax3,'Units', 'normalized', 'Position', [0.43,0.13,.2,.1]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.6 0.45]);
colormap(handles.cmap)

%% Opening and closing callbacks

function resetGUI_Callback(hObject, eventdata, handles)
% Loads new GUI session
Bayesian;

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Closes current GUI session
MODAclose(hObject,handles)
