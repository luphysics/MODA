%Version 1.01
%**************************************************************************
%***************************** CoherenceMulti GUI ******************************
%**************************************************************************

%---------------------------Credits----------------------------------------
% Wavelet Transform: Dmytro Iatsenko
% Wavelet phase coherence: Dmytro Iatsenko

%----------------------------Documentation---------------------------------
% Calculates the wavelet transforms of signal pairs, loaded as .csv or .mat
% files. Then calculates the wavelet phase coherence using phases from the
% wavelet transform.



function varargout = CoherenceMulti(varargin)
% COHERENCEMULTI MATLAB code for CoherenceMulti.fig
%      COHERENCEMULTI, by itself, creates a new COHERENCEMULTI or raises the existing
%      singleton*.
%
%      H = COHERENCEMULTI returns the handle to a new COHERENCEMULTI or the handle to
%      the existing singleton*.
%
%      COHERENCEMULTI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COHERENCEMULTI.M with the given input arguments.
%
%      COHERENCEMULTI('Property','Value',...) creates a new COHERENCEMULTI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CoherenceMulti_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CoherenceMulti_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help CoherenceMulti

% Last Modified by GUIDE v2.5 22-Mar-2018 14:55:03
%*************************************************************************%
%                BEGIN initialization code - DO NOT EDIT                  %
%                ----------------------------------------                 %
%*************************************************************************%
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CoherenceMulti_OpeningFcn, ...
                   'gui_OutputFcn',  @CoherenceMulti_OutputFcn, ...
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


function CoherenceMulti_OpeningFcn(hObject, eventdata, handles, varargin)
handles=MODAsettings(hObject, handles);

% Disable plot functions on startup
set(handles.plot_TS,'Enable','off')
set(handles.save_3dplot,'Enable','off')
set(handles.save_both_plot,'Enable','off')
set(handles.save_avg_plot,'Enable','off')
set(handles.save_mm_plot,'Enable','off')

% Disable save functions on startup
set(handles.save_avg_csv,'Enable','off')
set(handles.save_avg_mat,'Enable','off')

drawnow;
handles.output = hObject;
guidata(hObject, handles);


function file_read_Callback(hObject, eventdata, handles)
% Loads user data when selected from File --> Load time series

[handles,A]=MODAreadcheck(handles);
    if A==1
%------- Resetting axes and clearing data from previous ---
    %------- calculations. ------------------------------------
    cla(handles.plot3d,'reset');
    cla(handles.plot_pow,'reset');
    cla(handles.cum_avg,'reset');
    cla(handles.time_series_1,'reset');
    cla(handles.time_series_2,'reset');
    uistack(handles.plot_pow,'top');
    uistack(handles.plot3d,'top');
    set(handles.plot3d,'visible','on');
    set(handles.plot_pow,'visible','on');
    set(handles.cum_avg,'visible','off');
    
    if isfield(handles, 'freqarr');handles = rmfield(handles, 'freqarr');else end
    if isfield(handles, 'sig');handles = rmfield(handles, 'sig');else end
    if isfield(handles, 'time_avg_wpc');handles = rmfield(handles, 'time_avg_wpc');else end
    if isfield(handles, 'leg1');handles = rmfield(handles, 'leg1');else end
    if isfield(handles, 'time_axis');handles = rmfield(handles, 'time_axis');else end
    if isfield(handles, 'time_axis_ds');handles = rmfield(handles, 'time_axis_ds');else end
    if isfield(handles, 'TPC_surr_avg_arr');handles = rmfield(handles, 'TPC_surr_avg_arr');else end
    if isfield(handles, 'TPC_surr_avg_max');handles = rmfield(handles, 'TPC_surr_avg_max');else end
    if isfield(handles, 'TPC');handles = rmfield(handles, 'TPC');else end
    if isfield(handles, 'surrogates');handles = rmfield(handles, 'surrogates');else end
    if isfield(handles, 'leg');handles = rmfield(handles, 'leg');else end
    if isfield(handles, 'wopt');handles = rmfield(handles, 'wopt');else end
    if isfield(handles, 'sampling_freq');handles = rmfield(handles, 'sampling_freq');else end
    if isfield(handles, 'peak_value');handles = rmfield(handles, 'peak_value');else end
    if isfield(handles, 'thresh');handles = rmfield(handles, 'thresh');else end
    %--------------------------------------------
       
    
    linkaxes([handles.time_series_1 handles.time_series_2],'x');
    set(handles.signal_list,'Value',1);
    
    % Load data
    [handles,sig]=MODAread(handles,1,"even");
    
    if sig==0
    else
    % Create signal list
    list = cell(size(sig,1)/2+1,1);
    list{1,1} = 'Signal Pair 1';
    for i = 2:size(sig,1)/2
        list{i,1} = sprintf('Signal Pair %d',i);
    end
    set(handles.signal_list,'String',list);
    list{size(sig,1)/2+1,1} = sprintf('Average Plot (All)');
    set(handles.signal_list,'String',list);     
    
   
    % Plot time series
    plot(handles.time_series_1,handles.time_axis,handles.sig(1,:),'color',handles.linecol(1,:));
    xlim(handles.time_series_1,[0,size(sig,2)./handles.sampling_freq]);
    plot(handles.time_series_2,handles.time_axis,handles.sig(1+size(sig,1)/2,:),'color',handles.linecol(1,:));
    xlim(handles.time_series_2,[0,size(sig,2)./handles.sampling_freq]);    
    refresh_limits_Callback(hObject, eventdata, handles);
    guidata(hObject,handles);  
    ylabel(handles.time_series_1,'Sig 1');
    ylabel(handles.time_series_2,'Sig 2');
    xlabel(handles.time_series_2,'Time (s)');
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');
    set(handles.signal_length,'String',strcat(num2str(size(sig,2)/handles.sampling_freq/60),' minutes'));
    set(handles.time_series_1,'ytickmode','auto','yticklabelmode', 'auto','xticklabel',[],'FontUnits','normalized');
    set(handles.time_series_2,'ytickmode','auto','yticklabelmode', 'auto','FontUnits','normalized');
    end
    else
        return;
    end


function varargout = CoherenceMulti_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function status_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','Please Import Signal');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function status_Callback(hObject, eventdata, handles, msg)
set(handles.status,'String',msg);
drawnow;

function intervals_Callback(hObject, eventdata, handles)
% Marking lines on the graphs    
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
    if interval_selected == size(handles.sig,1)/2 + 1 
        xl = get(child_handles(i),'ylim');
        for j = 1:size(intervals,2)            
            x = [xl(1) xl(2)];        
            z = ones(1,size(x,2));
            y = intervals(j)*ones(1,size(x,2));
            plot3(handles.cum_avg,y,x,z,'--k');
            xticks = get(handles.cum_avg,'xtick');
            xticks = unique(sort([xticks intervals]));
            set(handles.cum_avg,'xtick',xticks);            
        end        
        set(child_handles(i),'ylim',xl);       
        
        
    elseif length(interval_selected)>1
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
        set(child_handles(i),'ylim',xl);
    else       
        if(size(intervals)>0)
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

            end    
        end
    end
    
    %set(handles.plot_pow,'Yticklabel',[]);

function wavlet_transform_Callback(hObject, eventdata, handles)
handles.currsig=[];
try

handles.st=get(handles.surrogate_type,'Value');
handles.nsurr=str2double(get(handles.surrogate_count,'String'));
handles=MODAwpc(hObject, eventdata, handles, 1);

delete(handles.h)
xyplot_Callback(hObject, eventdata, handles);
intervals_Callback(hObject, eventdata, handles)
catch e
    errordlg(e.message,'Error');
    set(handles.wt_single,'Enable','on')
    set(handles.wavlet_transform,'Enable','on')
    
    rethrow(e)
end
guidata(hObject,handles);
%end

% set(handles.wt_single,'Enable','off')
% set(handles.wavlet_transform,'Enable','off')
% set(handles.subtract_surrogates,'Enable','off')
% try
% 
%     fmax = str2double(get(handles.max_freq,'String'));
%     fs = handles.sampling_freq;
%     f0 =  str2double(get(handles.central_freq,'String')); A=f0<=0.4;
%     items = get(handles.wavelet_type,'String');
%     index_selected = get(handles.wavelet_type,'Value');
%     wtype = items{index_selected}; B=strcmp(wtype,'Bump');
% 
%     if (A+0)+(B+0)==2
%           errordlg('The bump wavelet requires that f0 > 0.4. Please enter a higher value.','Parameter Error');
%           set(handles.wt_single,'Enable','on')
%           set(handles.wavlet_transform,'Enable','on')
%           return;
%     end
%     
%     if fmax>fs/2
%           errordlg(['Maximum frequency cannot be higher than the Nyquist frequency. Please enter a value less than or equal to ',num2str(fs/2),' Hz.'],'Parameter Error');
%           set(handles.wt_single,'Enable','on')
%           set(handles.wavlet_transform,'Enable','on')
%           return;
%     end
%     
%     status_Callback(hObject, eventdata, handles, 'Calculating Wavelet Transform...');
%     fs = handles.sampling_freq;
%     fmin = str2double(get(handles.min_freq,'String'));
%     fc =  str2double(get(handles.central_freq,'String'));
%     surrogate_count = floor(str2double(get(handles.surrogate_count,'String')));
%     
%     if (surrogate_count == 1)
%         errordlg('Number of surrogates must be greater than 1','Parameter Error');
%         set(handles.status,'String','Enter Valid Parameters before continuing');
%         drawnow;
%         return;
%     end
%     
%     if isnan(fs)
%       errordlg('Sampling frequency must be specified','Parameter Error');
%       set(handles.status,'String','Enter Valid Parameters before continuing');
%       drawnow;
%       return;
%     end
%     
%     if ~isfield(handles,'sig')
%       errordlg('Signal not found','Signal Error');
%       set(handles.status,'String','Enter Valid Parameters before continuing');
%       drawnow;
%       return;
%     end
%     
%     items = get(handles.wavelet_type,'String');
%     index_selected = get(handles.wavelet_type,'Value');
%     wavelet_type_selected = items{index_selected};
%     
%     items = get(handles.preprocess,'String');
%     index_selected = get(handles.preprocess,'Value');
%     preprocess_selected = items{index_selected};
%     
%     items = get(handles.cutedges,'String');
%     index_selected = get(handles.cutedges,'Value');
%     cutedges_selected = items{index_selected};
%     
%     
%     sig = handles.sig;    
%     
%     xl = csv_to_mvar(get(handles.xlim,'String'));
%     xl = xl.*fs;
%     xl(2) = min(xl(2),size(sig,2));
%     xl(1) = max(xl(1),1);
%     xl = xl./fs;
%     time_axis = xl(1):1/fs:xl(2);
%     
%     if length(time_axis)>=2000
%         screensize = max(get(groot,'Screensize'));
%         under_sample = floor(size(sig,2)/screensize);
%     else
%         under_sample = 1;
%     end
% 
%     handles.time_axis_ds = time_axis(1:under_sample:end);
%     n = size(handles.sig,1)/2 ;
%     
% % Taking only selected part of the signal
%     xl = get(handles.xlim,'String');
%     xl = csv_to_mvar(xl);
%     xl = xl.*fs;
%     xl(2) = min(xl(2),size(handles.sig,2));
%     xl(1) = max(xl(1),1);
%     sig = sig(:,xl(1):xl(2));
%     
%    
%         if fmin<=1/(length(handles.sig)/fs)
%           errordlg(['WT minimum frequency too low. To automatically calculate for minimum possible frequency leave "Min Freq" field blank.'],'Parameter Error'); 
%           set(handles.wt_single,'Enable','on')
%           set(handles.wavlet_transform,'Enable','on')
%           return;
%         end
%     
%     
%     set(handles.status,'String','Calculating Wavelet Transform...');
%     
%     handles.h = waitbar(0,'Calculating coherence...',...
%             'CreateCancelBtn',...
%             'setappdata(gcbf,''canceling'',1)');
%     setappdata(handles.h,'canceling',0)
%     guidata(hObject,handles);
%     handles.surimp=1;
%    
%     for p = 1:n
%         if getappdata(handles.h,'canceling')
%             handles.surimp=0;
%             guidata(hObject,handles);
%             set(handles.wt_single,'Enable','on')
%             set(handles.wavlet_transform,'Enable','on')
%             break;
%         else
%         end
%         
%         status_Callback(hObject, eventdata, handles, sprintf('Calculating Wavelet Transform of Signal %d/%d',p,n));
%         if(isnan(fmax)&& isnan(fmin))
%             if(isnan(fc))                              
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);      
%             end
%         elseif(isnan(fmax))
%             if(isnan(fc))
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%             end
%         elseif(isnan(fmin))
%             if(isnan(fc))
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%             end
%         else
%             if(isnan(fc))
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%             end
%         end
%         
%         handles.TPC{p,1} = tlphcoh(wt_1,wt_2,handles.freqarr,fs);
%         %handles.time_avg_wpc{p,1} = nanmean(handles.TPC{p,1}.');
%         handles.time_avg_wpc{p,1} = wphcoh(wt_1,wt_2);
%         handles.TPC{p,1} = handles.TPC{p,1}(:,1:under_sample:end);
%         waitbar(p/n,handles.h);
%         
%         if getappdata(handles.h,'canceling')
%             handles.surimp=0;
%             guidata(hObject,handles);
%             set(handles.wt_single,'Enable','on')
%             set(handles.wavlet_transform,'Enable','on')
%             break;
%         end
%     end   
%     delete(handles.h);
%     %---------------
%             
%     % Surrogate Calculation
%  
%     items = get(handles.surrogate_type,'String');
%     index_selected = get(handles.surrogate_type,'Value');
%     surrogate_type = items{index_selected};
%     handles.stype=surrogate_type;
%     if (surrogate_count > 1) && (floor(surrogate_count) == surrogate_count)
%         handles.h = waitbar(0,'Calculating surrogates...',...
%             'CreateCancelBtn',...
%             'setappdata(gcbf,''canceling'',1)');
%         setappdata(handles.h,'canceling',0)
%         guidata(hObject,handles);
%         
%         if getappdata(handles.h,'canceling')
%                 set(handles.wt_single,'Enable','on')
%                 set(handles.wavlet_transform,'Enable','on')
%                 return;
%             
%         end
%         if handles.surimp==0;
%             set(handles.wt_single,'Enable','on')
%             set(handles.wavlet_transform,'Enable','on')
%             return;
%         end
%         handles.surrogates = cell(size(handles.sig,1),1);
%         se=size(handles.sig,1)/2;
%         for p = 1:se
%             status_Callback(hObject, eventdata, handles, ['Calculating surrogates for signal ',num2str(p),' of ',num2str(se)]); 
%             handles.surrogates{p,1} = surrcalc(sig(p+se,:),surrogate_count,surrogate_type,0,fs); 
%         end    
%         status_Callback(hObject, eventdata, handles, 'Surrogates complete');
%         handles.TPC_surr_avg_arr = cell(surrogate_count,size(handles.sig,1)/2);
%         
%         
% c=1;
%         for idx = 1:size(handles.sig,1)/2
%             if getappdata(handles.h,'canceling')
%                 set(handles.wt_single,'Enable','on')
%                 set(handles.wavlet_transform,'Enable','on')
%                 break;
%             
%             end
%             status_Callback(hObject, eventdata, handles, sprintf('Calculating Wavelet Transform of Signal %d/%d',idx,size(handles.sig,1)/2));
%             if(isnan(fmax)&& isnan(fmin))
%                 if(isnan(fc))                              
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);     
%                 end
%             elseif(isnan(fmax))
%                 if(isnan(fc))
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%                 end
%             elseif(isnan(fmin))
%                 if(isnan(fc))
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%                 end
%             else
%                 if(isnan(fc))
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);
%                 end
%             end
%             for p = 1:surrogate_count
%                 if getappdata(handles.h,'canceling')
%                     errordlg(['Calculation interrupted by user. Only ',num2str(p-1),' surrogate(s) calculated'],'Warning');
%                     set(handles.wt_single,'Enable','on')
%                     set(handles.wavlet_transform,'Enable','on')
%                     break;
%                 end
%                 status_Callback(hObject, eventdata, handles, ...
%                     sprintf('Calculating Wavelet Transform for surrogate:%d/%d for signal pair %d',p,surrogate_count,idx));
% 
%                 if(isnan(fmax)&& isnan(fmin))
%                     if(isnan(fc))                  
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);     
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);    
%                     end
%                 elseif(isnan(fmax))
%                     if(isnan(fc))
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);  
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
%                     end
%                 elseif(isnan(fmin))
%                     if(isnan(fc))
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);    
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
%                     end
%                 else
%                     if(isnan(fc))
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
%                     end
%                 end
% 
%                 TPC_surrogate = tlphcoh(wt_1,WT_surrogate,handles.freqarr,fs);
%                 %handles.TPC_surr_avg_arr{p,idx} = nanmean(TPC_surrogate.');
%                 handles.TPC_surr_avg_arr{p,idx} = wphcoh(wt_1,WT_surrogate);
%                 c=c+1;
%                 waitbar(c/(surrogate_count*(size(handles.sig,1)/2)),handles.h);
%                 if getappdata(handles.h,'canceling')
%                     errordlg(['Calculation interrupted by user. Only ',num2str(p),' surrogate(s) calculated'],'Warning');
%                     set(handles.wt_single,'Enable','on')
%                     set(handles.wavlet_transform,'Enable','on')
%                     break;
%                 end
%             end
%         delete(handles.h)    
%         end
%         
% 
%         status_Callback(hObject, eventdata, handles, 'Finished calculating surrogates');
%     %-------------------------------------------------------------------------      
% 
%         surrogate_analysis = get(handles.surrogate_analysis,'Value');        
%         handles.TPC_surr_avg_max = cell(size(handles.sig,1)/2,1);
%         alph = str2double(get(handles.surrogate_percentile,'String'));
%         handles.thresh=surrogate_analysis;
%         for i = 1:size(handles.sig,1)/2
%             if(surrogate_analysis == 2)
%                 t = cell2mat(handles.TPC_surr_avg_arr);
%                 t = t(:,length(handles.freqarr)*(i-1)+1:length(handles.freqarr)*(i));
%                 
%                 if floor((surrogate_count+1)*alph)==0
%                     handles.TPC_surr_avg_max{i,1} = max(t);                
%                 else
%                     K=floor((surrogate_count+1)*alph);
%                     s1=sort(t,'descend');
%                     handles.TPC_surr_avg_max{i,1}= s1(K,:);              
%                
%                 end
%                 
%  
%             elseif(surrogate_analysis == 1)    
%                 t = cell2mat(handles.TPC_surr_avg_arr);
%                 t = t(:,length(handles.freqarr)*(i-1)+1:length(handles.freqarr)*(i));
%                 handles.TPC_surr_avg_max{i,1} = max(t);
% 
%             end
%         end     
%     end
%     guidata(hObject,handles);
%     xyplot_Callback(hObject, eventdata, handles);
%     intervals_Callback(hObject, eventdata, handles)
%     guidata(hObject,handles);
%     set(handles.intervals,'Enable','on')
%     set(handles.wt_single,'Enable','on')
%     set(handles.plot_TS,'Enable','on')
%     set(handles.save_3dplot,'Enable','on')
%     set(handles.save_both_plot,'Enable','on')
%     set(handles.save_avg_plot,'Enable','on')
%     set(handles.save_mm_plot,'Enable','on')
%     set(handles.save_avg_csv,'Enable','on')
%     set(handles.save_avg_mat,'Enable','on')
%     set(handles.wavlet_transform,'Enable','on')
%     if surrogate_count > 1
%     set(handles.subtract_surrogates,'Enable','on')
%     else
%     end
%     drawnow;
%     
% catch e
%     errordlg(e.message,'Error');
%     set(handles.wt_single,'Enable','on')
%     set(handles.wavlet_transform,'Enable','on')
%     if surrogate_count > 1
%     set(handles.subtract_surrogates,'Enable','on')
%     else
%     end
%     delete(handles.h)
%     rethrow(e)
% end

function wt_single_Callback(hObject, eventdata, handles)
handles.currsig=get(handles.signal_list,'Value');
try
handles=MODAwpc(hObject, eventdata, handles, 0);
delete(handles.h)
xyplot_Callback(hObject, eventdata, handles);
intervals_Callback(hObject, eventdata, handles)
catch e
    errordlg(e.message,'Error');
    set(handles.wt_single,'Enable','on')
    set(handles.wavlet_transform,'Enable','on')
    
    rethrow(e)
end
guidata(hObject,handles);
% handles.time_avg_wpc=[];
% set(handles.wt_single,'Enable','off')
% set(handles.wavlet_transform,'Enable','off')
% set(handles.subtract_surrogates,'Enable','off')
% try
%     
%     fmax = str2double(get(handles.max_freq,'String'));
%     fs = handles.sampling_freq;
%     f0 =  str2double(get(handles.central_freq,'String')); A=f0<=0.4;
%     items = get(handles.wavelet_type,'String');
%     index_selected = get(handles.wavelet_type,'Value');
%     wtype = items{index_selected}; B=strcmp(wtype,'Bump');
% 
%     if (A+0)+(B+0)==2
%           errordlg('The bump wavelet requires that f0 > 0.4. Please enter a higher value.','Parameter Error');
%           set(handles.wt_single,'Enable','on')
%           set(handles.wavlet_transform,'Enable','on')
%           return;
%     end
%     
%     if fmax>fs/2
%           errordlg(['Maximum frequency cannot be higher than the Nyquist frequency. Please enter a value less than or equal to ',num2str(fs/2),' Hz.'],'Parameter Error');
%           set(handles.wt_single,'Enable','on')
%           set(handles.wavlet_transform,'Enable','on')
%           return;
%     end
% 
%     status_Callback(hObject, eventdata, handles, 'Calculating Wavelet Transform...');
%     fs = handles.sampling_freq;
%     fmin = str2double(get(handles.min_freq,'String'));
%     fc =  str2double(get(handles.central_freq,'String'));
%     surrogate_count = floor(str2double(get(handles.surrogate_count,'String')));
%     
%     if (surrogate_count == 1)
%         errordlg('Number of surrogates must be greater than 1','Parameter Error');
%         set(handles.status,'String','Enter Valid Parameters before continuing');
%         drawnow;
%         return;
%     end
%     
%     if isnan(fs)
%       errordlg('Sampling frequency must be specified','Parameter Error');
%       set(handles.status,'String','Enter Valid Parameters before continuing');
%       drawnow;
%       return;
%     end
%     
%     if ~isfield(handles,'sig')
%       errordlg('Signal not found','Signal Error');
%       set(handles.status,'String','Enter Valid Parameters before continuing');
%       drawnow;
%       return;
%     end
%     
%     items = get(handles.wavelet_type,'String');
%     index_selected = get(handles.wavelet_type,'Value');
%     wavelet_type_selected = items{index_selected};
%     
%     items = get(handles.preprocess,'String');
%     index_selected = get(handles.preprocess,'Value');
%     preprocess_selected = items{index_selected};
%     
%     items = get(handles.cutedges,'String');
%     index_selected = get(handles.cutedges,'Value');
%     cutedges_selected = items{index_selected};
%     
%     
%     sig = handles.sig;    
%     
%     xl = csv_to_mvar(get(handles.xlim,'String'));
%     xl = xl.*fs;
%     xl(2) = min(xl(2),size(sig,2));
%     xl(1) = max(xl(1),1);
%     xl = xl./fs;
%     time_axis = xl(1):1/fs:xl(2);
%     
%     if length(time_axis)>=2000
%         screensize = max(get(groot,'Screensize'));
%         under_sample = floor(size(sig,2)/screensize);
%     else 
%         under_sample = 1;
%     end
% 
%     handles.time_axis_ds = time_axis(1:under_sample:end);
%     n = size(handles.sig,1)/2 ;
%     
% % Taking only selected part of the signal
%     xl = get(handles.xlim,'String');
%     xl = csv_to_mvar(xl);
%     xl = xl.*fs;
%     xl(2) = min(xl(2),size(handles.sig,2));
%     xl(1) = max(xl(1),1);
%     sig = sig(:,xl(1):xl(2));
%     
%     if fmin<=1/(length(handles.sig)/fs)
%           errordlg(['WT minimum frequency too low. To automatically calculate for minimum possible frequency leave "Min Freq" field blank.'],'Parameter Error'); 
%           set(handles.wt_single,'Enable','on')
%           set(handles.wavlet_transform,'Enable','on')
%           return;
%     end
%     
%     set(handles.status,'String','Calculating Wavelet Transform...');
%     
%      handles.h = waitbar(0,'Calculating coherence...',...
%             'CreateCancelBtn',...
%             'setappdata(gcbf,''canceling'',1)');
%     setappdata(handles.h,'canceling',0)
%     guidata(hObject,handles);
%     handles.surimp=1;
%     for p = get(handles.signal_list,'Value')
%         
%         
%         if getappdata(handles.h,'canceling')
%             handles.surimp=0;
%             guidata(hObject,handles);
%             set(handles.wt_single,'Enable','on')
%             set(handles.wavlet_transform,'Enable','on')
%             break;
%         end
%         
%         status_Callback(hObject, eventdata, handles, sprintf('Calculating Wavelet Transform of Signal %d/%d',p,n));
%         if(isnan(fmax)&& isnan(fmin))
%             if(isnan(fc))                              
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);      
%             end
%         elseif(isnan(fmax))
%             if(isnan(fc))
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%             end
%         elseif(isnan(fmin))
%             if(isnan(fc))
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%             end
%         else
%             if(isnan(fc))
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%             else
%                     [wt_1,handles.freqarr,handles.wopt]=wt(sig(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
% 
%                     [wt_2,handles.freqarr,handles.wopt]=wt(sig(p+n,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                         'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%             end
%         end
%         
%         handles.TPC{p,1} = tlphcoh(wt_1,wt_2,handles.freqarr,fs);
%         %handles.time_avg_wpc{p,1} = nanmean(handles.TPC{p,1}.');
%         handles.time_avg_wpc{p,1} = wphcoh(wt_1,wt_2);
%         handles.TPC{p,1} = handles.TPC{p,1}(:,1:under_sample:end);
%         waitbar(p/p,handles.h);
%         
%         if getappdata(handles.h,'canceling')
%             handles.surimp=0;
%             guidata(hObject,handles);
%             set(handles.wt_single,'Enable','on')
%             set(handles.wavlet_transform,'Enable','on')
%             break;
%         end
%         
%     end   
%     delete(handles.h)
%     %---------------
%     
%     % Surrogate Calculation
%     
%     
%      
%     items = get(handles.surrogate_type,'String');
%     index_selected = get(handles.surrogate_type,'Value');
%     surrogate_type = items{index_selected};
%     handles.stype=surrogate_type;
%     if (surrogate_count > 1) && (floor(surrogate_count) == surrogate_count)
%         handles.h = waitbar(0,'Calculating surrogates...',...
%             'CreateCancelBtn',...
%             'setappdata(gcbf,''canceling'',1)');
%         setappdata(handles.h,'canceling',0)
%         guidata(hObject,handles);
%         if handles.surimp==0;
%             set(handles.wt_single,'Enable','on')
%             set(handles.wavlet_transform,'Enable','on')
%             return;
%         end
%         handles.surrogates = cell(size(handles.sig,1),1);
%         se=size(handles.sig,1)/2;
%         for p =get(handles.signal_list,'Value')
%             status_Callback(hObject, eventdata, handles, ['Calculating surrogates for signal ',num2str(p),' of ',num2str(se)]);
%             handles.surrogates{p,1} = surrcalc(sig(p+se,:),surrogate_count,surrogate_type,0,fs);   
%             status_Callback(hObject, eventdata, handles, 'Surrogates complete');
%         end    
%                 
%         
% 
%         for idx = get(handles.signal_list,'Value')
%             if getappdata(handles.h,'canceling')                
%                 set(handles.wt_single,'Enable','on')
%                 set(handles.wavlet_transform,'Enable','on')
%                 break;
%             end
%             status_Callback(hObject, eventdata, handles, sprintf('Calculating Wavelet Transform of Signal %d/%d',idx,size(handles.sig,1)/2));
%             if(isnan(fmax)&& isnan(fmin))
%                 if(isnan(fc))                              
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);     
%                 end
%             elseif(isnan(fmax))
%                 if(isnan(fc))
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%                 end
%             elseif(isnan(fmin))
%                 if(isnan(fc))
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc); 
%                 end
%             else
%                 if(isnan(fc))
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected); 
%                 else
%                         [wt_1,handles.freqarr]=wt(sig(idx,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);
%                 end
%             end
%             for p = 1:surrogate_count
%                 if getappdata(handles.h,'canceling')
%                     errordlg(['Calculation interrupted by user. Only ',num2str(p-1),' surrogate(s) calculated'],'Warning');
%                     set(handles.wt_single,'Enable','on')
%                     set(handles.wavlet_transform,'Enable','on')
%                     break;
%                 end
%                 status_Callback(hObject, eventdata, handles, ...
%                     sprintf('Calculating Wavelet Transform for surrogate:%d/%d for signal pair %d',p,surrogate_count,idx));
% 
%                 if(isnan(fmax)&& isnan(fmin))
%                     if(isnan(fc))                  
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);     
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);    
%                     end
%                 elseif(isnan(fmax))
%                     if(isnan(fc))
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);  
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
%                     end
%                 elseif(isnan(fmin))
%                     if(isnan(fc))
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);    
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
%                     end
%                 else
%                     if(isnan(fc))
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected);
%                     else
%                         [WT_surrogate,handles.freqarr]=wt(handles.surrogates{idx,1}(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutedges_selected,...
%                             'Preprocess',preprocess_selected,'Wavelet',wavelet_type_selected,'f0',fc);  
%                     end
%                 end
% 
%                 TPC_surrogate = tlphcoh(wt_1,WT_surrogate,handles.freqarr,fs);
%                 %handles.TPC_surr_avg_arr{p,idx} = nanmean(TPC_surrogate.'); 
%                 handles.TPC_surr_avg_arr{p,idx} = wphcoh(wt_1,WT_surrogate); 
%                 waitbar(p/surrogate_count,handles.h);
%                 if getappdata(handles.h,'canceling')
%                     errordlg(['Calculation interrupted by user. Only ',num2str(p),' surrogate(s) calculated'],'Warning');
%                     guidata(hObject,handles);
%                     set(handles.wt_single,'Enable','on')
%                     set(handles.wavlet_transform,'Enable','on')
%                 break;
%                 end
%             end
%         end
% delete(handles.h)
%         status_Callback(hObject, eventdata, handles, 'Finished calculating surrogates');
%     %-------------------------------------------------------------------------      
% 
%         surrogate_analysis = get(handles.surrogate_analysis,'Value');        
%         handles.TPC_surr_avg_max = cell(size(handles.sig,1)/2,1);
%         handles.thresh=surrogate_analysis;
%         for i = get(handles.signal_list,'Value'); 
%             if(surrogate_analysis == 2)
%                 t = cell2mat(handles.TPC_surr_avg_arr);
%                 t = t(:,1:length(handles.freqarr));        
%                 surrogate_percentile = str2double(get(handles.surrogate_percentile,'String'));
%                 handles.TPC_surr_avg_max{i,1} = prctile(t,surrogate_percentile);
% 
%             elseif(surrogate_analysis == 1)    
%                 t = cell2mat(handles.TPC_surr_avg_arr);
%                 t = t(:,1:length(handles.freqarr));
%                 handles.TPC_surr_avg_max{i} = max(t);
% 
%             end
%         end     
%     end
%     guidata(hObject,handles);
%     xyplot_Callback(hObject, eventdata, handles);
%     intervals_Callback(hObject, eventdata, handles)
%     guidata(hObject,handles);
%     set(handles.intervals,'Enable','on')
%     set(handles.wt_single,'Enable','on')
%     set(handles.plot_TS,'Enable','on')
%     set(handles.save_3dplot,'Enable','on')
%     set(handles.save_both_plot,'Enable','on')
%     set(handles.save_avg_plot,'Enable','on')
%     set(handles.save_mm_plot,'Enable','on')
%     set(handles.save_avg_csv,'Enable','on')
%     set(handles.save_avg_mat,'Enable','on')
%     set(handles.wavlet_transform,'Enable','on')
%     if surrogate_count > 1
%     set(handles.subtract_surrogates,'Enable','on')
%     else
%     end
%     drawnow;
% catch e
%     errordlg(e.message,'Error');
%     set(handles.wt_single,'Enable','on')
%     set(handles.wavlet_transform,'Enable','on')
%     if surrogate_count > 1
%     set(handles.subtract_surrogates,'Enable','on')
%     else
%     end
%     delete(handles.h)
%     rethrow(e)
% end


function xyplot_Callback(hObject, eventdata, handles)
% Plots all figures
    signal_selected = get(handles.signal_list,'Value');  
    if isfield(handles,'time_avg_wpc') && signal_selected~=(size(handles.sig_cut,1)/2)+1 && isempty(handles.time_avg_wpc{signal_selected,1})
        cla(handles.plot3d,'reset')
        cla(handles.plot_pow,'reset')
    else
    
    if any(signal_selected == size(handles.sig,1)/2+1) && isfield(handles,'freqarr')     
        cla(handles.plot3d,'reset');
        cla(handles.plot_pow,'reset');
        cla(handles.cum_avg,'reset');
        set(handles.plot3d,'visible','off');
        set(handles.plot_pow,'visible','off');
        set(handles.cum_avg,'visible','on');
        set(handles.save_3dplot,'Enable','off')
        set(handles.save_both_plot,'Enable','off')
        set(handles.save_avg_plot,'Enable','off')
        set(handles.save_mm_plot,'Enable','on')
        hold(handles.cum_avg,'on');
        uistack(handles.cum_avg, 'top');
        
        if size(handles.sig,1)/2 > 1
            plot(handles.cum_avg, handles.freqarr, mean(cell2mat(handles.time_avg_wpc)),'-','Linewidth',3,'color',handles.linecol(1,:));
            plot(handles.cum_avg, handles.freqarr, median(cell2mat(handles.time_avg_wpc)),'--','Linewidth',3,'color',handles.linecol(2,:));
        else
            plot(handles.cum_avg, handles.freqarr, cell2mat(handles.time_avg_wpc),'-','Linewidth',3,'color',handles.linecol(1,:));
            plot(handles.cum_avg, handles.freqarr, cell2mat(handles.time_avg_wpc),'--','Linewidth',3,'color',handles.linecol(2,:));
        end
        ylabel(handles.cum_avg,'Overall Coherence','FontUnits','Points','Fontsize',10);
        xlabel(handles.cum_avg,'Frequency (Hz)','FontUnits','Points','Fontsize',10);
        handles.leg1={'Mean','Median'};
        legend(handles.cum_avg,handles.leg1)
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
            if(signal_selected(i) <= size(handles.sig,1)/2)                                
                plot(handles.cum_avg, handles.freqarr, handles.time_avg_wpc{signal_selected(i),1},sty,'color',handles.linecol(ind,:),'LineWidth',handles.line2width);                     
                ylabel(handles.cum_avg,'Overall Coherence','FontUnits','Points','Fontsize',10);
                xlabel(handles.cum_avg,'Frequency (Hz)','FontUnits','Points','Fontsize',10);            
                [M,I] = max(handles.time_avg_wpc{signal_selected(i),1});
                handles.leg1{i+2}=['Pair ',num2str(signal_selected(i))];
                legend(handles.cum_avg,handles.leg1)

            end
        end
        
        set(handles.cum_avg,'xscale','log');     
        idx_first = find(sum(~isnan(handles.time_avg_wpc{1,1}),1) > 0, 1 ,'first');
        idx_last = find(sum(~isnan(handles.time_avg_wpc{1,1}),1) > 0, 1 , 'last');   
        xlim(handles.cum_avg,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
        grid(handles.cum_avg,'off');
        box(handles.cum_avg,'on');
        
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
        
        handles.peak_value = max(handles.TPC{signal_selected,1}(:))+.1;
        pcolor(handles.plot3d, handles.time_axis_ds , handles.freqarr, handles.TPC{signal_selected,1});    
        
        plot(handles.plot_pow, handles.time_avg_wpc{signal_selected,1}, handles.freqarr,'LineWidth',2,'color',handles.linecol(1,:));
        set(handles.plot3d,'FontUnits','Points','Fontsize',10);
        set(handles.plot_pow,'FontUnits','Points','Fontsize',10);
        hold(handles.plot_pow,'on');
        if str2double(get(handles.surrogate_count,'String'))>1            
            plot(handles.plot_pow,handles.TPC_surr_avg_max{signal_selected,1} , handles.freqarr,'LineWidth',2,'color',handles.linecol(2,:));
            handles.leg = {'Original Signal','Surrogate'};
            pow_plot_leg = legend(handles.plot_pow, handles.leg);
            set(pow_plot_leg,'fontsize',8,'Position',[0.855 0.92 0.05 0.05]);
        else
        end
        hold(handles.plot_pow,'off');
        xlabel(handles.plot_pow,'Overall Coherence','FontUnits','Points','Fontsize',10);       
        
        
        xlabel(handles.plot3d,'Time (s)','FontUnits','Points','Fontsize',10);
        ylabel(handles.plot3d,'Frequency (Hz)','FontUnits','Points','Fontsize',10);    
        ylabel(handles.plot_pow,'Frequency (Hz)','FontUnits','Points','Fontsize',10);
        set(handles.plot3d,'FontUnits','Points','Fontsize',10);
        set(handles.plot_pow,'FontUnits','Points','Fontsize',10);
        
        colormap(handles.plot3d,handles.cmap);
        shading(handles.plot3d,'interp');       
        set(handles.plot3d,'yscale','log');
        set(handles.plot_pow,'yscale','log');%,'yticklabel',[]);        
        set(handles.plot3d,'ylim',[min(handles.freqarr) max(handles.freqarr)],...
            'xlim',[handles.time_axis_ds(1) handles.time_axis_ds(end)]);
        
        idx_first = find(sum(~isnan(handles.time_avg_wpc{signal_selected,1}),1) > 0, 1 ,'first');
        idx_last = find(sum(~isnan(handles.time_avg_wpc{signal_selected,1}),1) > 0, 1 , 'last');      
        ylim(handles.plot_pow,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
        ylim(handles.plot3d,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
        set(handles.status,'String','Done Plotting'); 
        grid(handles.plot3d,'on');
        grid(handles.plot_pow,'off');
        
    end
    
    set(handles.plot3d,'FontUnits','Points','Fontsize',10);
    set(handles.plot_pow,'FontUnits','Points','Fontsize',10);
    set(handles.cum_avg,'FontUnits','Points','Fontsize',10);
    
    guidata(hObject,handles);
    end

%---------------------------Surrogate Analysis-----------------
function handles=surrogate_analysis_Callback(hObject, eventdata, handles)

surrogate_analysis = get(handles.surrogate_analysis,'Value');  
if(surrogate_analysis == 2)
    set(handles.surrogate_percentile,'Enable','on');
elseif(surrogate_analysis == 1)
    set(handles.surrogate_percentile,'Enable','off');
else
end

function surrplot(hObject, eventdata, handles)
surrogate_analysis = get(handles.surrogate_analysis,'Value');        
handles.TPC_surr_avg_max = cell(size(handles.sig_cut,1)/2,1);

        
        alph = str2double(get(handles.surrogate_percentile,'String'));
        handles.thresh=surrogate_analysis;
        
        if length(cell2mat(handles.TPC_surr_avg_arr))>length(handles.freqarr)
        for i = 1:size(handles.sig_cut,1)/2
            if(surrogate_analysis == 2)
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,length(handles.freqarr)*(i-1)+1:length(handles.freqarr)*(i));
                
                if floor((handles.nscalc+1)*alph)==0
                    handles.TPC_surr_avg_max{i,1} = max(t);                
                else 
                    K=floor((handles.nscalc+1)*alph);
                    s1=sort(t,'descend');
                    handles.TPC_surr_avg_max{i,1}= s1(K,:);              
               
                end
                
 
            elseif(surrogate_analysis == 1)   && handles.nscalc>1 
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,length(handles.freqarr)*(i-1)+1:length(handles.freqarr)*(i));
                handles.TPC_surr_avg_max{i,1} = max(t);

            end
        end   
        else
            for i = get(handles.signal_list,'Value'); 
            if(surrogate_analysis == 2)
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,1:length(handles.freqarr));   
                
                if floor((handles.nscalc+1)*alph)==0
                    handles.TPC_surr_avg_max{i,1} = max(t);                
                else 
                    K=floor((handles.nscalc+1)*alph);
                    s1=sort(t,'descend');
                    handles.TPC_surr_avg_max{i,1}= s1(K,:);              
               
                end
                

            elseif(surrogate_analysis == 1)    
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,1:length(handles.freqarr));
                handles.TPC_surr_avg_max{i} = max(t);

            end
            end     
            
        end


guidata(hObject,handles);
xyplot_Callback(hObject, eventdata, handles);
drawnow;

subtract_surrogates_Callback(hObject, eventdata, handles)
guidata(hObject,handles);


% function surrogate_percentile_Callback(hObject, eventdata, handles)
% set(handles.surrogate_percentile,'Enable','on');
% handles=surrogate_analysis_Callback(hObject, eventdata, handles);
% guidata(hObject,handles);
% subtract_surrogates_Callback(hObject, eventdata, handles)
% guidata(hObject,handles);


function subtract_surrogates_Callback(hObject, eventdata, handles)
    if str2double(get(handles.surrogate_count,'String'))<2
    else
    
signal_selected = get(handles.signal_list,'Value');
toggle = get(handles.subtract_surrogates,'Value');
if toggle == get(handles.subtract_surrogates,'Max')
    cla(handles.plot_pow);
    corrected_coherence = handles.time_avg_wpc{signal_selected,1} - handles.TPC_surr_avg_max{signal_selected,1};
    corrected_coherence = subplus(corrected_coherence);
    plot(handles.plot_pow ,corrected_coherence, handles.freqarr,'LineWidth',2,'color',handles.linecol(1,:));
    set(handles.plot_pow,'yscale','log');%,'yticklabel',[]);     
    idx_first = find(sum(~isnan(handles.time_avg_wpc{signal_selected,1}),1) > 0, 1 ,'first');
    idx_last = find(sum(~isnan(handles.time_avg_wpc{signal_selected,1}),1) > 0, 1 , 'last');      
    ylim(handles.plot_pow,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
    xlabel(handles.plot_pow,'Overall Coherence','FontUnits','Points','Fontsize',10);
    ylabel(handles.plot_pow,'Frequency (Hz)','FontUnits','Points','Fontsize',10); 
    handles.leg = {'Surrogate Subtracted'};
    pow_plot_leg = legend(handles.plot_pow, handles.leg);
    set(pow_plot_leg,'fontsize',8,'Position',[0.855 0.92 0.05 0.05]);
    
else
    cla(handles.plot_pow);
    hold(handles.plot_pow,'on');
    plot(handles.plot_pow ,handles.time_avg_wpc{signal_selected,1}, handles.freqarr,'LineWidth',2,'color',handles.linecol(1,:));
    if(size(handles.TPC_surr_avg_max)>0)
        plot(handles.plot_pow ,handles.TPC_surr_avg_max{signal_selected,1} , handles.freqarr,'LineWidth',2,'color',handles.linecol(2,:));
    end
    hold(handles.plot_pow,'off');     
    set(handles.plot_pow,'yscale','log');%,'yticklabel',[]);     
    idx_first = find(sum(~isnan(handles.time_avg_wpc{signal_selected,1}),1) > 0, 1 ,'first');
    idx_last = find(sum(~isnan(handles.time_avg_wpc{signal_selected,1}),1) > 0, 1 , 'last');      
    ylim(handles.plot_pow,[handles.freqarr(idx_first) handles.freqarr(idx_last)]);
    set(handles.status,'String','Done Plotting');
    xlabel(handles.plot_pow,'Overall Coherence','FontUnits','Points','Fontsize',10);   
    ylabel(handles.plot_pow,'Frequency (Hz)','FontUnits','Points','Fontsize',10); 
    handles.leg = {'Original Signal','Surrogate'};
    pow_plot_leg = legend(handles.plot_pow, handles.leg);
    set(pow_plot_leg,'fontsize',8,'Position',[0.855 0.92 0.05 0.05]);
    
end
    end
grid(handles.plot_pow,'off');
intervals_Callback(hObject, eventdata, handles)
set(handles.plot_pow,'FontUnits','Points','Fontsize',10);
guidata(hObject,handles);

function detrend_signal_popup_Callback(hObject, eventdata, handles)
cla(handles.plot_pp,'reset');
    

function signal_list_Callback(hObject, eventdata, handles)

    signal_selected = get(handles.signal_list, 'Value');
    
    if any(signal_selected == size(handles.sig,1)/2+1)
        set(handles.signal_list,'Max',size(handles.sig,1)/2);
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
    
    if any(signal_selected ~= size(handles.sig,1)/2+1) && length(signal_selected) == 1
        
        plot(handles.time_series_1, handles.time_axis, handles.sig(signal_selected,:),'color',handles.linecol(1,:));
        xl = csv_to_mvar(get(handles.xlim, 'String'));
        xlim(handles.time_series_1, xl);
                
        plot(handles.time_series_2, handles.time_axis, handles.sig(signal_selected+size(handles.sig,1)/2,:),'color',handles.linecol(1,:));
        xlim(handles.time_series_2, xl);        
        
        refresh_limits_Callback(hObject, eventdata, handles);
        set(handles.status, 'String', 'Select Data And Continue With Wavelet Transform');
        if isfield(handles,'TPC')
            xyplot_Callback(hObject, eventdata, handles);
        end
        intervals_Callback(hObject, eventdata, handles)
        ylabel(handles.time_series_1,'Sig 1'); 
        ylabel(handles.time_series_2,'Sig 2'); 
        xlabel(handles.time_series_2, 'Time (s)');
        set(handles.time_series_1,'xticklabel',[],'ytickmode','auto','yticklabelmode', 'auto','FontUnits','normalized');
        set(handles.time_series_2,'ytickmode','auto','yticklabelmode', 'auto','FontUnits','normalized');
    elseif any(signal_selected == size(handles.sig,1)/2+1)
        xyplot_Callback(hObject, eventdata, handles);
        intervals_Callback(hObject, eventdata, handles)
    end

% --------------------------------------------------------------------
    
%---------------------------Limits-----------------------------
function xlim_Callback(hObject, eventdata, handles)
% When the values of xlim are changed the graphs are updated
    xl = csv_to_mvar(get(handles.xlim,'String'));
    xlim(handles.time_series_1,xl);
    xlim(handles.time_series_2,xl);
    xlim(handles.plot_pp,xl);
    t = xl(2) - xl(1);
    set(handles.length,'String',t);

function ylim_Callback(hObject, eventdata, handles)
% When the values of ylim are changed the graphs are updated  
    yl = csv_to_mvar(get(handles.ylim,'String'));
    ylim(handles.time_series_1,yl);
    ylim(handles.time_series_2,yl);


%---------------------------Updating Value of limits Limits-----------------------------
function refresh_limits_Callback(hObject, eventdata, handles)
% Calculates limits of the plot    
    
    x = get(handles.time_series_1,'xlim');
    t = x(2) - x(1);
    x = strcat([num2str(x(1)),', ',num2str(x(2))]);    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat([num2str(y(1)),', ',num2str(y(2))]);
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);
    
    

% ---------------------------Zoom Updating--------------------------
function zoom_in_OffCallback(hObject, eventdata, handles)
% Refreshes the limit values right after the tool is deselected
    x = get(handles.time_series_1,'xlim');
    t = x(2) - x(1);
    x = strcat([num2str(x(1)),', ',num2str(x(2))]);    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat([num2str(y(1)),', ',num2str(y(2))]);
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);

% -----------------------------Zoom Updating--------------------------
function zoom_out_OffCallback(hObject, eventdata, handles)
% Refreshes the limit values right after the tool is deselected
    x = get(handles.time_series_1,'xlim');
    t = x(2) - x(1);
    x = strcat([num2str(x(1)),', ',num2str(x(2))]);    
    
    y = get(handles.time_series_1,'ylim');
    y = strcat([num2str(y(1)),', ',num2str(y(2))]);
    
    set(handles.xlim,'String',x);
    set(handles.ylim,'String',y);
    set(handles.length,'String',t);
    
function plot_type_SelectionChangeFcn(hObject, eventdata, handles)

    switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
        case 'power'
            plot_type = 1;
        case 'amp'
            plot_type = 2;
    end
    data = guidata(hObject);
    data.plot_type = plot_type;
    guidata(hObject,data); 


function plot_TS_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.time_series_1, Fig);
ax2 = copyobj(handles.time_series_2, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.1,0.55,.85,.35],'FontUnits','points','FontSize',10);
set(ax2,'Units', 'normalized', 'Position', [0.1,0.15,.85,.35],'FontUnits','points','FontSize',10);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);


function save_3dplot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.plot3d, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
colormap(ax,handles.cmap); 
h = colorbar;
ylabel(h, 'Wavelet coherence')


function save_avg_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.plot_pow, Fig);
view(90,-90);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
legend(ax,handles.leg,'location','best')

function save_both_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax1 = copyobj(handles.plot3d, Fig);
h = colorbar;
ylabel(h, 'Wavelet coherence')
colormap(Fig,handles.cmap);
ax2 = copyobj(handles.plot_pow, Fig);
set(ax1,'Units', 'normalized', 'Position', [0.07,0.2,.55,.7]);
set(ax2,'Units', 'normalized', 'Position', [0.8,0.2,.18,.7], 'YTickMode', 'auto', 'YTickLabelMode', 'auto');
%ylabel(ax2,[])
set(Fig,'Units','normalized','Position', [0.2 0.2 0.6 0.5]);


function save_mm_plot_Callback(hObject, eventdata, handles)
Fig = figure;
ax = copyobj(handles.cum_avg, Fig);
set(ax,'Units', 'normalized', 'Position', [0.1,0.2,.85,.7]);
set(Fig,'Units','normalized','Position', [0.2 0.2 0.5 0.5]);
legend(handles.leg1)

function save_avg_csv_Callback(hObject, eventdata, handles)
try
[FileName,PathName] = uiputfile('.csv','Save Coherence Data');
if FileName==0
    return;
else
end
save_location = strcat(PathName,FileName);
avg_coh = cell2mat(handles.time_avg_wpc)';
xl = csv_to_mvar(get(handles.xlim,'String'));
L=xl(2)*handles.wopt.fs - xl(1)*handles.wopt.fs;

Coherence_data.Coherence=avg_coh;
Coherence_data.Frequency=handles.freqarr;
if str2double(get(handles.surrogate_count,'String'))>1
Coherence_data.Surrogates=cell2mat(handles.TPC_surr_avg_max)';
else
end
Coherence_data.Time=linspace(xl(1),xl(2),L);
Coherence_data.Sampling_frequency=handles.wopt.fs;
Coherence_data.fmax=handles.wopt.fmax;
Coherence_data.fmin=handles.wopt.fmin;
Coherence_data.fr=handles.wopt.f0;
Coherence_data.Preprocessing=handles.wopt.Preprocess;
Coherence_data.Cut_Edges=handles.wopt.CutEdges;
Coherence_data.Wavelet_type=handles.wopt.Wavelet;

if str2double(get(handles.surrogate_count,'String'))>1
Coherence_data.Surrogate_type=handles.stype;
Coherence_data.Surrogate_number=get(handles.surrogate_count,'String');
if handles.thresh==2
    Coherence_data.Surrogate_threshold=['Significance ', get(handles.surrogate_percentile,'String')];
else
    Coherence_data.Surrogate_threshold='Maximum';
end
%     if strcmp(get(handles.surrogate_analysis,'String'),'Percentile')
%         Coherence_data.Surrogate_type=get(handles.surrogate_percentile,'String');
%     else
%     end
else
end
data=csvsaving(Coherence_data,handles);
cell2csv(save_location,data,',');
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end


function save_avg_mat_Callback(hObject, eventdata, handles)
try
[FileName,PathName] = uiputfile('.mat','Save Coherence Data');
if FileName==0
    return;
else
end
save_location = strcat(PathName,FileName)
avg_coh = cell2mat(handles.time_avg_wpc)';
xl = csv_to_mvar(get(handles.xlim,'String'));
L=xl(2)*handles.wopt.fs - xl(1)*handles.wopt.fs;

Coherence_data.Coherence=avg_coh;
if ~isempty(handles.currsig)
    Coherence_data.Selected_sig=handles.currsig;
else
end
Coherence_data.Frequency=handles.freqarr;
if str2double(get(handles.surrogate_count,'String'))>1
Coherence_data.Surrogates=cell2mat(handles.TPC_surr_avg_max)';
else
end
Coherence_data.Time=linspace(xl(1),xl(2),L);
Coherence_data.Sampling_frequency=handles.wopt.fs;
Coherence_data.fmax=handles.wopt.fmax;
Coherence_data.fmin=handles.wopt.fmin;
Coherence_data.fr=handles.wopt.f0;
Coherence_data.Preprocessing=handles.wopt.Preprocess;
Coherence_data.Cut_Edges=handles.wopt.CutEdges;
Coherence_data.Wavelet_type=handles.wopt.Wavelet;

if str2double(get(handles.surrogate_count,'String'))>1
Coherence_data.Surrogate_type=handles.stype;
Coherence_data.Surrogate_number=get(handles.surrogate_count,'String');
if handles.thresh==2
    Coherence_data.Surrogate_threshold=['Significance ', get(handles.surrogate_percentile,'String')];
else
    Coherence_data.Surrogate_threshold='Maximum';
end
%     if strcmp(get(handles.surrogate_analysis,'String'),'Percentile')
%         Coherence_data.Surrogate_type=get(handles.surrogate_percentile,'String');
%     else
%     end
else
end

save(save_location,'Coherence_data');
catch e
    errordlg(e.message,'Error')
    rethrow(e)
end

function data=csvsaving(D,handles)
L=length(D.Frequency);
N=size(D.Coherence,2);

if isfield(D,'Surrogates')
    data=cell(L+13,(N*2)+1);
    dstart=18;
    data{1,1}='MODA v1.0 - Wavelet Phase Coherence';
    data{2,1}=date;
    data{3,1}=[];
    data{4,1}='PARAMETERS';
    data{5,1}='Sampling frequency (Hz)';
    data{5,2}=D.Sampling_frequency;
    data{6,1}='Maximum frequency (Hz)';
    data{6,2}=D.fmax;
    data{7,1}='Minimum frequency (Hz)';
    data{7,2}=D.fmin;
    data{8,1}='Central frequency';
    data{8,2}=D.fr;
    data{9,1}='Preprocessing';
    data{9,2}=D.Preprocessing;
    data{10,1}='Wavelet type';
    data{10,2}=D.Wavelet_type;
    data{11,1}='Cut Edges';
    data{11,2}=D.Cut_Edges;
    data{12,1}='Time start (s)';
    data{12,2}=min(D.Time);
    data{13,1}='Time end (s)';
    data{13,2}=max(D.Time);
    data{14,1}='Surrogate type';
    data{14,2}=D.Surrogate_type;
    data{15,1}='Surrogate number';
    data{15,2}=D.Surrogate_number;
    data{16,1}='Surrogate threshold';
    data{16,2}=D.Surrogate_threshold;
else
    data=cell(L+13,N+1);
    dstart=15;
    data{1,1}='Wavelet phase coherence toolbox';
    data{2,1}=date;
    data{3,1}=[];
    data{4,1}='PARAMETERS';
    data{5,1}='Sampling frequency (Hz)';
    data{5,2}=D.Sampling_frequency;
    data{6,1}='Maximum frequency (Hz)';
    data{6,2}=D.fmax;
    data{7,1}='Minimum frequency (Hz)';
    data{7,2}=D.fmin;
    data{8,1}='Central frequency';
    data{8,2}=D.fr;
    data{9,1}='Preprocessing';
    data{9,2}=D.Preprocessing;
    data{10,1}='Wavelet type';
    data{10,2}=D.Wavelet_type;
    data{11,1}='Cut Edges';
    data{11,2}=D.Cut_Edges;
    data{12,1}='Time start (s)';
    data{12,2}=min(D.Time);
    data{13,1}='Time end (s)';
    data{13,2}=max(D.Time);
    
end


data{dstart,1}='Frequency';
for l=1:L;
data{l+dstart,1}=D.Frequency(l);
end

if isfield(D,'Surrogates')

for j=1:N
    if isempty(handles.currsig)
        data{dstart,j+1}=['Coherence ',num2str(j)];
        data{dstart,j+N+1}=['Surrogate ',num2str(j)];
    else
        data{dstart,j+1}=['Coherence ',num2str(handles.currsig)];
        data{dstart,j+N+1}=['Surrogate ',num2str(handles.currsig)];
    end
    for k=1:L
        data{k+dstart,j+1}=D.Coherence(k,j);  
        data{k+dstart,j+N+1}=D.Surrogates(k,j);
    end   
    
end

else
    
for j=1:N
    if isempty(handles.currsig)
        data{dstart,j+1}=['Coherence ',num2str(j)];
    else
        data{dstart,j+1}=['Coherence ',num2str(handles.currsig)];
    end
      for k=1:L
        data{k+dstart,j+1}=D.Coherence(k,j);  
      end   
    
end
end


% --------------------------------------------------------------------


% --------------------------------------------------------------------
function resetGUI_Callback(hObject, eventdata, handles)

CoherenceMulti;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

MODAclose(hObject,handles)


% --------------------------------------------------------------------
function save_session_Callback(hObject, eventdata, handles)

MODAsave(handles)


% --------------------------------------------------------------------
function load_session_Callback(hObject, eventdata, handles)

handles=MODAload;


% --- Executes on button press in supdate.
function supdate_Callback(hObject, eventdata, handles)

stype=get(handles.surrogate_type,'Value');

nsurr=str2double(get(handles.surrogate_count,'String'));
save stype stype nsurr

if isfield(handles,'nsurr') && stype==handles.st && nsurr==handles.nsurr
    surrplot(hObject, eventdata, handles)
elseif isfield(handles,'nsurr') && stype==handles.st && nsurr~=handles.nsurr
    errordlg('Number of surrogates changed, please recalculate coherence','Error')
elseif isfield(handles,'nsurr') && stype~=handles.st && nsurr==handles.nsurr
    errordlg('Surrogate type changed, please recalculate coherence','Error')
end
