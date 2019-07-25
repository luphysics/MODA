% Wavelet phase coherence calculation in MODA

function handles=MODAwpc(hObject, eventdata, handles, type)
% type - 1=multiple, 2= single
try
handles.time_avg_wpc=cell(size(handles.sig_cut,1),1);
handles.nscalc=(str2double(get(handles.surrogate_count,'String')));
set(handles.wt_single,'Enable','off')
set(handles.wavlet_transform,'Enable','off')
set(handles.subtract_surrogates,'Enable','off')

set(handles.status,'String', 'Calculating Wavelet phase coherence..');

fmax = str2double(get(handles.max_freq,'String'));
fmin = str2double(get(handles.min_freq,'String'));
f0 =  str2double(get(handles.central_freq,'String')); A=f0<=0.4;
fs = handles.sampling_freq;
fc =  str2double(get(handles.central_freq,'String'));
items = get(handles.wavelet_type,'String'); index_selected = get(handles.wavelet_type,'Value'); wtype = items{index_selected}; B=strcmp(wtype,'Bump');
ns=(str2double(get(handles.surrogate_count,'String')));
handles.TPC_surr_avg_arr=cell(size(handles.sig_cut,1),ns);
items = get(handles.surrogate_type,'String'); index_selected = get(handles.surrogate_type,'Value');
surrogate_type = items{index_selected}; handles.stype=surrogate_type;

if index_selected==1
    surrogate_type='RP';
else
end



    if (A+0)+(B+0)==2
          errordlg('The bump wavelet requires that f0 > 0.4. Please enter a higher value.','Parameter Error');
          set(handles.wt_single,'Enable','on')
          set(handles.wavlet_transform,'Enable','on')
          return;
    end
    
    if fmax>fs/2
          errordlg(['Maximum frequency cannot be higher than the Nyquist frequency. Please enter a value less than or equal to ',num2str(fs/2),' Hz.'],'Parameter Error');
          set(handles.wt_single,'Enable','on')
          set(handles.wavlet_transform,'Enable','on')
          return;
    end    
    
    if (ns == 1)
        errordlg('Number of surrogates must be greater than 1','Parameter Error');
        set(handles.status,'String','Enter Valid Parameters before continuing');
        drawnow;
        return;
    end
    
    if isnan(fs)
      errordlg('Sampling frequency must be specified','Parameter Error');
      set(handles.status,'String','Enter Valid Parameters before continuing');
      drawnow;
      return;
    end
    
    if ~isfield(handles,'sig')
      errordlg('Signal not found','Signal Error');
      set(handles.status,'String','Enter Valid Parameters before continuing');
      drawnow;
      return;
    end
    
    if fmin<=1/(length(handles.sig)/fs)
          errordlg(['WT minimum frequency too low. To automatically calculate for minimum possible frequency leave "Min Freq" field blank.'],'Parameter Error'); 
          set(handles.wt_single,'Enable','on')
          set(handles.wavlet_transform,'Enable','on')
          return;
    end
    
    
    items = get(handles.preprocess,'String'); index_selected = get(handles.preprocess,'Value');  ppselect = items{index_selected};
    
    items = get(handles.cutedges,'String'); index_selected = get(handles.cutedges,'Value'); cutselect = items{index_selected};
    
    if length(handles.sig_cut)>=2000
        screensize = max(get(groot,'Screensize'));
        under_sample = floor(size(handles.sig_cut,2)/screensize);
    else 
        under_sample = 1;
    end

    handles.time_axis_ds = handles.time_axis_cut(1:under_sample:end);
    n = size(handles.sig_cut,1)/2 ;
    
    handles.h = waitbar(0,'Calculating coherence...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
    setappdata(handles.h,'canceling',0)
    guidata(hObject,handles);
        
    % Choose whether to calculate for all signals or only current signal
    if type==1
        inds=1:n;
    else
        inds=get(handles.signal_list,'Value');
    end
    
    handles.surrogates = cell(size(handles.sig_cut,1)/2,1);
    handles.TPC_surr_avg_arr = cell(ns,size(handles.sig_cut,1)/2);
    
    % Calculate coherence (and surrogates if required)
    completed=0;
    for p = inds
        if ishandle(handles.h)
        if  getappdata(handles.h,'canceling')
            set(handles.wt_single,'Enable','on')
            set(handles.wavlet_transform,'Enable','on')
            delete(handles.h)
            guidata(hObject,handles);
            break;
        else
        end
        else
            break;
        end
        
        set(handles.status,'String', sprintf('Calculating wavelet phase coherence of Signal %d/%d',p,n));
        [wt_1,handles]=wtcalc(handles.sig_cut(p,:),handles,fmax,fmin,fc,ppselect,cutselect,wtype,fs);
        [wt_2,handles]=wtcalc(handles.sig_cut(p+n,:),handles,fmax,fmin,fc,ppselect,cutselect,wtype,fs);
        
        handles.TPC{p,1} = tlphcoh(wt_1,wt_2,handles.freqarr,fs);
        handles.time_avg_wpc{p,1} = wphcoh(wt_1,wt_2);
        handles.TPC{p,1} = handles.TPC{p,1}(:,1:under_sample:end);
        
        if (ns > 1) %&& (floor(ns) == ns)
            if ishandle(handles.h)
            if getappdata(handles.h,'canceling')
                set(handles.wt_single,'Enable','on')
                set(handles.wavlet_transform,'Enable','on')
                delete(handles.h)
                guidata(hObject,handles);
                return;
            else
            end
            else
                break;
            end
            
            set(handles.status,'String', ['Calculating surrogates for signal ',num2str(p),' of ',num2str(n)]);
            
            %if strcmp(ppselect,'off')
                handles.surrogates{p,1} = surrcalc(handles.sig_cut(p+n,:),ns,surrogate_type,0,fs); 
                           
              
            for k=1:ns
                pause(0.00001)
                if ishandle(handles.h)
                if getappdata(handles.h,'canceling')
                    set(handles.wt_single,'Enable','on')
                    set(handles.wavlet_transform,'Enable','on')
                    delete(handles.h)
                    guidata(hObject,handles);
                    break;
                else
                end
                else
                    break;
                end
                [WT_surrogate,handles]=wtcalc(handles.surrogates{p,1}(k,:),handles,fmax,fmin,fc,ppselect,cutselect,wtype,fs);
                
                handles.TPC_surr_avg_arr{k,p} = wphcoh(wt_1,WT_surrogate);
                
            end
           
        end
        
       
        
        
        if ishandle(handles.h)
            
        waitbar(p/n,handles.h);
        else
            
        end
        if p==inds(end)
            completed=1;   
        else
        end
    end   
    guidata(hObject, handles);
    
    delete(handles.h)
    
    if completed==1
    
        surrogate_analysis = get(handles.surrogate_analysis,'Value');        
        handles.TPC_surr_avg_max = cell(size(handles.sig_cut,1)/2,1);
        
        alph = str2double(get(handles.surrogate_percentile,'String'));
        handles.thresh=surrogate_analysis;
        
        if type==1
        for i = 1:size(handles.sig_cut,1)/2
            if(surrogate_analysis == 2)
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,length(handles.freqarr)*(i-1)+1:length(handles.freqarr)*(i));
                
                if floor((ns+1)*alph)==0
                    handles.TPC_surr_avg_max{i,1} = max(t);                
                else 
                    K=floor((ns+1)*alph);
                    s1=sort(t,'descend');
                    handles.TPC_surr_avg_max{i,1}= s1(K,:);              
               
                end
                
 
            elseif(surrogate_analysis == 1)   && ns>1 
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
                
                if floor((ns+1)*alph)==0
                    handles.TPC_surr_avg_max{i,1} = max(t);                
                else 
                    K=floor((ns+1)*alph);
                    s1=sort(t,'descend');
                    handles.TPC_surr_avg_max{i,1}= s1(K,:);              
               
                end
                

            elseif(surrogate_analysis == 1)  && ns>1 
                t = cell2mat(handles.TPC_surr_avg_arr);
                t = t(:,1:length(handles.freqarr));
                handles.TPC_surr_avg_max{i} = max(t);

            end
            end     
            
        end
        
    guidata(hObject, handles);
        
    set(handles.intervals,'Enable','on')
    set(handles.wt_single,'Enable','on')
    set(handles.plot_TS,'Enable','on')
    set(handles.save_3dplot,'Enable','on')
    set(handles.save_both_plot,'Enable','on')
    set(handles.save_avg_plot,'Enable','on')
    set(handles.save_mm_plot,'Enable','on')
    set(handles.save_avg_csv,'Enable','on')
    set(handles.save_avg_mat,'Enable','on')
    set(handles.wavlet_transform,'Enable','on')
    if ns > 1
        set(handles.subtract_surrogates,'Enable','on')
    else
    end
    
    set(handles.status,'String', sprintf('Calculation complete'))
    else
        set(handles.status,'String', sprintf('Calculation interrupted by user'))
        set(handles.wt_single,'Enable','on')
        set(handles.wavlet_transform,'Enable','on')
    end
catch e
    errordlg(e.message,'Error');
    set(handles.wt_single,'Enable','on')
    set(handles.wavlet_transform,'Enable','on')
    
    rethrow(e)
end
    
    function [wt_1,handles]=wtcalc(sig,handles,fmax,fmin,fc,ppselect,cutselect,wtype,fs)
        
        
        if(isnan(fmax)&& isnan(fmin))
            if(isnan(fc))                              
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype); 

               
            else
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc); 

                        
            end
        elseif(isnan(fmax))
            if(isnan(fc))
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'fmin',fmin,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype); 

                    
            else
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'fmin',fmin,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc); 

                   
            end
        elseif(isnan(fmin))
            if(isnan(fc))
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype); 

                   
            else
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc); 

                   
            end
        else
            if(isnan(fc))
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype); 

                  
            else
                    [wt_1,handles.freqarr,handles.wopt]=wt(sig,fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc); 

                   
            end
        end

    

