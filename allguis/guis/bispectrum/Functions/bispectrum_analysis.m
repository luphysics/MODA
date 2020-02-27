function [amp1, amp2, avg_amp1, avg_amp2, pow1, pow2, avg_pow1, avg_pow2, bispxxx, bispppp, bispxpp, bisppxx, surrxxx, surrppp, surrxpp, surrpxx] = bispectrum_analysis(sig1, sig2, fmax, fmin, fc, fs, nv, ns, preprocess)

bispxxx=[];
bispppp=[];
bispxpp=[];
bisppxx=[];
surrxxx=[];
surrppp=[];
surrxpp=[];
surrpxx=[];

if ~preprocess
    [bispxxx] = bispecWavNew(sig1,sig1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',1,'wbar',1); bispxxx=abs(bispxxx);
    [bispppp] = bispecWavNew(sig2,sig2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',2,'wbar',1); bispppp=abs(bispppp);
    [bispxpp,freqarr, wavopt,WT1, WT2] = bispecWavNew(sig1,sig2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',3,'wbar',1); bispxpp=abs(bispxpp);
    [bisppxx] = bispecWavNew(sig2,sig1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',4,'wbar',1); bisppxx=abs(bisppxx);
    
    if ns>0
        for j=1:ns
            surr1=wavsurrogate(sig1,'IAAFT2',1);
            surr2=wavsurrogate(sig2,'IAAFT2',1);
            
            surrxxx(:,:,j)=abs(bispecWavNew(surr1,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',1));
            surrppp(:,:,j)=abs(bispecWavNew(surr2,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',2));
            surrxpp(:,:,j)=abs(bispecWavNew(surr1,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',3));
            surrpxx(:,:,j)=abs(bispecWavNew(surr2,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'preprocess','off','handles',handles,'hObject',hObject,'num',4));
        end
    end
else
    [bispxxx] = bispecWavNew(sig1,sig1,fs,'handles',handles,'hObject',hObject,'f0',fc,'nv',nv,'fmin',fmin,'fmax',fmax,'num',1,'wbar',1); bispxxx=abs(bispxxx);
    [bispppp] = bispecWavNew(sig2,sig2,fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'handles',handles,'hObject',hObject,'num',2,'wbar',1); bispppp=abs(bispppp);
    [bispxpp,freqarr,wavopt,WT1, WT2] = bispecWavNew(sig1,sig2,fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'handles',handles,'hObject',hObject,'num',3,'wbar',1); bispxpp=abs(bispxpp);
    [bisppxx] = bispecWavNew(sig2,sig1,fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'handles',handles,'hObject',hObject,'num',4,'wbar',1); bisppxx=abs(bisppxx);
    
    if ns>0
        for j=1:ns
            surr1=wavsurrogate(sig1,'IAAFT2',1);
            surr2=wavsurrogate(sig2,'IAAFT2',1);
            
            surrxxx(:,:,j)=abs(bispecWavNew(surr1,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'handles',handles,'hObject',hObject,'num',1));
            surrppp(:,:,j)=abs(bispecWavNew(surr2,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'handles',handles,'hObject',hObject,'num',2));
            surrxpp(:,:,j)=abs(bispecWavNew(surr1,surr2, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'handles',handles,'hObject',hObject,'num',3));
            surrpxx(:,:,j)=abs(bispecWavNew(surr2,surr1, fs,'fmin',fmin,'fmax',fmax,'f0',fc,'nv',nv,'handles',handles,'hObject',hObject,'num',4));
        end
    end
end

amp1=abs(WT1);
avg_amp1=nanmean(amp1,2);
pow1=abs(WT1).^2;
avg_pow1=nanmean(pow1,2);

amp2=abs(WT2);
avg_amp2=nanmean(amp2,2);
pow2=abs(WT2).^2;
avg_pow2=nanmean(pow2,2);

end

