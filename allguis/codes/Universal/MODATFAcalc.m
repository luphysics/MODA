function handles=MODATFAcalc(hObject, eventdata, handles,ty)
%ty - type of calculation (1=multiple, 2=single)
set(handles.plot_TS,'Enable','on')
set(handles.save_3dplot,'Enable','on')
set(handles.save_both_plot,'Enable','on')
set(handles.save_avg_plot,'Enable','on')
set(handles.save_mm_plot,'Enable','on')

set(handles.wt_single,'Enable','off')
set(handles.wavlet_transform,'Enable','off')
if ty==2
    set(handles.save_WT_coeff,'Enable','on')
elseif ty==1
    set(handles.save_WT_coeff,'Enable','off')
end
set(handles.save_session,'Enable','on')

handles.failed=false;
try
    % Obtain parameters from GUI
    fmax = str2double(get(handles.max_freq,'String'));
    fmin = str2double(get(handles.min_freq,'String'));
    fs = handles.sampling_freq;
    f0 =  str2double(get(handles.central_freq,'String')); A=f0<=0.4;
    items = get(handles.wavelet_type,'String');
    index_selected = get(handles.wavelet_type,'Value');
    wtype = items{index_selected}; B=strcmp(wtype,'Bump');    
    handles.fc =  str2double(get(handles.central_freq,'String'));

    if (A+0)+(B+0)==2
          errordlg('The bump wavelet requires that f0 > 0.4. Please enter a higher value.','Parameter Error');
          handles.failed = true;
          set(handles.wt_single,'Enable','on')
          set(handles.wavlet_transform,'Enable','on')
          return;
    end
    
    if fmax>fs/2
          errordlg(['Maximum frequency cannot be higher than the Nyquist frequency. Please enter a value less than or equal to ',num2str(fs/2),' Hz.'],'Parameter Error');
          handles.failed = true;
          set(handles.wt_single,'Enable','on')
          set(handles.wavlet_transform,'Enable','on')
          return;
    end
    
    %% Forces user to input minimum frequency for WFT, and changes resolution parameter according to fr/fmin, where fr is the user input resolution    
    if handles.calc_type==2 && isnan(fmin)
          errordlg(['Minimum frequency must be specified for WFT'],'Parameter Error');
          handles.failed = true;
          set(handles.wt_single,'Enable','on')
          set(handles.wavlet_transform,'Enable','on')
          return;
    elseif handles.calc_type==2
        handles.fc=f0/fmin;
    end
    
    if isnan(fs)
      errordlg('Sampling frequency must be specified','Parameter Error');
      handles.failed = true;
    end
    
    if handles.calc_type == 1
        set(handles.status,'String','Calculating Wavelet Transform...');
        %status_Callback(hObject, eventdata, handles, 'Calculating Wavelet Transform...');
    else
        set(handles.status,'String','Calculating Windowed Fourier Transform...');
        %status_Callback(hObject, eventdata, handles, 'Calculating Windowed Fourier Transform...');
    end 
    
    
    items = get(handles.wavelet_type,'String');
    index_selected = get(handles.wavelet_type,'Value');
    wtype = items{index_selected};
    
    if strcmp(wtype,'Kaiser')
        a=str2double(get(handles.kaisera,'String'));
        wtype = ['kaiser-',num2str(a)];        
    else
    end
    
    items = get(handles.preprocess,'String');
    index_selected = get(handles.preprocess,'Value');
    ppselect = items{index_selected};
    
    items = get(handles.cutedges,'String');
    index_selected = get(handles.cutedges,'Value');
    cutselect = items{index_selected};
    
    if ~isfield(handles,'sig')
      errordlg('Signal not found','Signal Error');
      handles.failed = true;
    end
    sig = handles.sig;    
    
    xl = csv_to_mvar(get(handles.xlim,'String'));
    xl = xl.*fs;
    xl(2) = min(xl(2),size(sig,2));
    xl(1) = max(xl(1),1);
    xl = xl./fs;
    time_axis = xl(1):1/fs:xl(2);
    if length(time_axis)>=2000
        screensize = max(get(groot,'Screensize'));
        under_sample = floor(size(sig,2)/screensize);
    else 
        under_sample = 1;
    end
    if handles.calc_type == 2
        under_sample = ceil(under_sample*3.5);
    end
    handles.time_axis_us = time_axis(1:under_sample:end);
    n = size(handles.sig,1) ;
    handles.WT = cell(n, 1);
    
% Taking only selected part of the signal
    xl = get(handles.xlim,'String');
    xl = csv_to_mvar(xl);
    xl = xl.*fs;
    xl(2) = min(xl(2),size(handles.sig,2));
    xl(1) = max(xl(1),1);
    handles.sig_cut = sig(:,xl(1):xl(2));
    
    
    if handles.calc_type==1
        if fmin<=1/(length(handles.sig_cut)/fs)
          errordlg(['WT minimum frequency too low. To automatically calculate for minimum possible frequency leave "Min Freq" field blank.'],'Parameter Error'); 
          handles.failed = true;
          set(handles.wt_single,'Enable','on')
          set(handles.wavlet_transform,'Enable','on')
          return;
        end
    else
    end
    if handles.calc_type == 1
        set(handles.status,'String','Calculating Wavelet Transform...');
    else
        set(handles.status,'String','Calculating Windowed Fourier Transform...');
    end
    
    handles.amp_WT = cell(n,1);
    handles.pow_WT = cell(n,1);
    handles.pow_arr = cell(n,1);
    handles.amp_arr = cell(n,1);
    
    handles.h = waitbar(0,'Calculating transform...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(handles.h,'canceling',0)
    guidata(hObject,handles);
    
    if ty==2
        I=get(handles.signal_list,'Value');
    else
        I=1:n;
    end
    
    for p = I
        if getappdata(handles.h,'canceling')
            break;
        end
        
        % Number of signals being transformed.
        if ty==1, count=n; 
        else, count=1; end
        
        if handles.calc_type == 1
            set(handles.status,'String',sprintf('Calculating Wavelet Transform of Signal %d/%d',p,count));
            
        else
            set(handles.status,'String',sprintf('Calculating Windowed Fourier Transform of Signal %d/%d',p,count));
            
        end
        wtwrapper;
        handles.WT=WT;
        WTamp = abs(WT);
        WTpow = abs(WT).^2;
        handles.pow_arr{p,1} = nanmean(WTpow.'); % Calculating Average Power
        handles.amp_arr{p,1} = nanmean(WTamp.'); % Calculating Average Amplitude  

        handles.amp_WT{p,1} = WTamp(:,1:under_sample:end);   
        handles.pow_WT{p,1} = WTpow(:,1:under_sample:end);
        waitbar(p/n,handles.h);
    end
    guidata(hObject,handles);
    
    delete(handles.h);
    set(handles.wt_single,'Enable','on')
    set(handles.wavlet_transform,'Enable','on')
    set(handles.mat_save,'Enable','on')
    set(handles.csv_save,'Enable','on')
    set(handles.save_WT_coeff,'Enable','on')
    set(handles.save_session,'Enable','on')
catch e
    errordlg(e.message,'Error');
    handles.failed = true;
    set(handles.wt_single,'Enable','on')
    set(handles.wavlet_transform,'Enable','on')
    delete(handles.h)
    rethrow(e)
end
guidata(hObject,handles);