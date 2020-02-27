
fs=handles.sampling_freq;
fc=handles.fc;
if(isnan(fmax)&& isnan(fmin))
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'Padding',0); 
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'Padding',0); 
                    end            
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc,'Padding',0); 
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'f0',fc,'Padding',0); 
                    end   
            end
        elseif(isnan(fmax))
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'Padding',0); 
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'Padding',0); 
                    end
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc,'Padding',0); 
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'f0',fc,'Padding',0); 
                    end
            end
        elseif(isnan(fmin))
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'Padding',0); 
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'Padding',0); 
                    end
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc,'Padding',0); 
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'f0',fc,'Padding',0); 
                    end
            end
        else
            if(isnan(fc))
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'Padding',0);
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'Padding',0); 
                    end
            else
                    if handles.calc_type == 1
                        [WT,handles.freqarr,handles.wopt]=wt(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Wavelet',wtype,'f0',fc,'Padding',0);
                    else
                        [WT,handles.freqarr,handles.wopt]=wft(handles.sig_cut(p,:),fs,'fmin',fmin,'fmax',fmax,'CutEdges',cutselect,...
                        'Preprocess',ppselect,'Window',wtype,'f0',fc,'Padding',0); 
                    end                 
            end
 end