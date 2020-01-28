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
                
                if(~isfield(handles, "fc"))
                   msg = "An error has occurred. Please re-calculate the transform before proceeding.";
                   errordlg(msg);
                   error(msg);
                end
                
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
                
                % This section relates to using "ridge" or "direct" reconstruction.
                %
                % Comments from Dymtro Iatsenko, paraphrased:
                %
                % [As far as I recall], the direct method for frequency is inferred directly from exact component reconstruction.
                %
                % Thus, in case of signal being a(t)cos(phi(t)) with slowly varying amplitude and frequency,
                % it will give exact amplitude, phase and frequency of this component, while the ridge method will not.
                %
                % Overall, the direct method is superior for signals where components are relatively clear;
                % however, with lots of noise and interferences the ridge method is more robust.
                
                % TO SWITCH FROM "direct" TO "ridge" RECONSTRUCTION, CHANGE 'direct' TO 'ridge' ON THE NEXT LINE.
                [handles.bands_iamp{j,k},handles.bands_iphi{j,k},handles.bands_freq{j,k}] = rectfr(tfsupp,WT,freqarr,wopt,'direct');
                
                handles.recon{j,k} = handles.bands_iamp{j,k}.*cos(handles.bands_iphi{j,k});
                handles.bands_iphi{j,k} = mod(handles.bands_iphi{j,k},2*pi);
                
                % Check if negative frequencies are present in the results.
                for i=1:size(handles.bands_freq, 1)
                    freq = handles.bands_freq{i};
                    negative = find(freq < 0);
                    
                    if ~isempty(negative)
                        msg = "Error: Negative frequencies are present in the result. This is an artefact of the direct reconstruction algorithm." + ...
                            newline + newline + ...
                            "To address this issue, please try one or more of the following:" + newline + ...
                            "1) Switching to the 'Lognorm' wavelet." + newline + ...
                            "2) Using a narrower frequency interval." + newline + ...
                            "3) Using a higher frequency resolution." + ...
                            newline + newline + ...
                            "For more information, please check the User Guide.";
                        
                        title = "Error";
                        see_docs = "Open User Guide";
                        dismiss = "Dismiss";
                        
                        % Show question dialog.
                        answer = questdlg(...
                            msg, ...
                            title, ...
                            see_docs, ...
                            dismiss, ...
                            dismiss ...
                            );
                        
                        if answer == see_docs
                            % Launch web browser with link to the relevant
                            % section of the User Guide.
                            web("https://github.com/luphysics/MODA/blob/master/docs/user-guide.md#negative-frequencies");
                        end
                        break;
                    end
                end
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