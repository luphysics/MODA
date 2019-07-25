function handles=MODAridge_filter(hObject,eventdata,handles)

list=get(handles.interval_list,'String');
if isempty(list)
    errordlg('Interval list is empty. Please select frequency bands for filtering','Error');
    return;
end
set(handles.transform,'Enable','off')
set(handles.filter_signal,'Enable','off')
set(handles.ridgecalc,'Enable','off')

%% Set up waitbar
    handles.h = waitbar(0,'Filtering...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(handles.h,'canceling',0)
    guidata(hObject,handles);

try
    list=get(handles.interval_list,'String');
    
    extype=handles.etype;
    
    if extype==2 %% If Butterworth filter is selected
        for j=1:size(handles.sig_cut,1)
            for k=1:size(list,1)
                fl=csv_to_mvar(list{k,1});
                [handles.bands{j,k},~] = loop_butter(handles.sig_cut(j,:),fl,handles.sampling_freq);
                handles.extract_phase{j,k}=angle(hilbert(handles.bands{j,k}));
                handles.extract_amp{j,k}=abs(hilbert(handles.bands{j,k}));
            end
        end
    elseif extype==1 %% If ridge extraction is selected
    % Window type
    wtypes=get(handles.wind_type,'String');
    wselect=get(handles.wind_type,'Value');
    wtype=wtypes{wselect};
    
     % Preprocessing input
    x=get(handles.preprocess,'String'); ind=get(handles.preprocess,'Value'); ppselect=x{ind}; 
    
    % Cut edges input
    x=get(handles.cutedges,'String'); ind=get(handles.cutedges,'Value'); cutselect=x{ind};
    for j=1:size(handles.sig_cut,1)
        for k=1:size(list,1)
            fl=csv_to_mvar(list{k,1});
           
    
    if(isnan(handles.fc))
         if handles.calc_type == 1
               [WT,freqarr,wopt]=wt(handles.sig_cut(j,:),handles.sampling_freq,'fmin',fl(1),'fmax',fl(2),'CutEdges','off',...
                'Preprocess',ppselect,'Wavelet',wtype);
         else
               [WT,freqarr,wopt]=wft(handles.sig_cut(j,:),handles.sampling_freq,'fmin',fl(1),'fmax',fl(2),'CutEdges','off',...
                'Preprocess',ppselect,'Window',wtype);
         end
    else
         if handles.calc_type == 1
               [WT,freqarr,wopt]=wt(handles.sig_cut(j,:),handles.sampling_freq,'fmin',fl(1),'fmax',fl(2),'CutEdges','off',...
                'Preprocess',ppselect,'Wavelet',wtype,'f0',handles.fc); 
         else
               [WT,freqarr,wopt]=wft(handles.sig_cut(j,:),handles.sampling_freq,'fmin',fl(1),'fmax',fl(2),'CutEdges','off',...
                'Preprocess',ppselect,'Window',wtype,'f0',handles.fc); 
         end
    end
                %Pre allocate for the cell structures
                tfsupp = ecurve(WT,freqarr,wopt);
                [handles.bands_iamp{j,k},handles.bands_iphi{j,k},handles.bands_freq{j,k}] = rectfr(tfsupp,WT,freqarr,wopt,'ridge');            
                handles.recon{j,k} = handles.bands_iamp{j,k}.*cos(handles.bands_iphi{j,k});
                handles.bands_iphi{j,k} = mod(handles.bands_iphi{j,k},2*pi);
        end
        waitbar(j/size(handles.sig_cut,1),handles.h)
    end
    
    
    else
        
    end
    set(handles.display_type,'Enable','on');
    set(handles.display_type,'Value',2);
    delete(handles.h);
    drawnow;
    guidata(hObject,handles);
    
    set(handles.transform,'Enable','on')
    set(handles.filter_signal,'Enable','on')
    set(handles.ridgecalc,'Enable','on')
    
    set(handles.save_filtered_sig_plot,'Enable','on')
    set(handles.save_ridge_plot,'Enable','on') 
    set(handles.save_phase_plot,'Enable','on') 
    set(handles.All_filt_plot,'Enable','on') 
    set(handles.save_csv,'Enable','on') 
    set(handles.save_mat,'Enable','on') 
    set(handles.save_session,'Enable','on')
    
    
catch e
    errordlg(e.message,'Error')
    set(handles.transform,'Enable','on')
    set(handles.filter_signal,'Enable','on')
    set(handles.ridgecalc,'Enable','on')
    delete(handles.h)
    rethrow(e)
    
end