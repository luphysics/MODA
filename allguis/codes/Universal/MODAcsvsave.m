function data=MODAcsvsave(D,ind)
%D=Filtered_data;
% ind=1 - filtered data
% ind=2 - filtered phases/freq
% ind=3 - filtered amplitudes

L=length(D.Time);
Nb=size(D.Freq_bands,1);

data{1,1}='MODA v1.0 - Filtering & Ridge Extraction';
data{2,1}=date;
data{3,1}=[];
data{4,1}='PARAMETERS';
data{5,1}='Sampling frequency (Hz)';
data{5,2}=D.Sampling_frequency;
data{6,1}='Filter type';
data{6,2}=D.Filter_type;
if isfield(D,'Preprocessing')
data{7,1}='Preprocessing';
data{7,2}=D.Preprocessing;
else
end
% data{8,1}='Cut Edges';
% data{8,2}=D.Cut_Edges;
data{9,1}='Time start (s)';
data{9,2}=min(D.Time);
data{10,1}='Time end (s)';
data{10,2}=max(D.Time);
for n=1:Nb
        data{n+10,1}=['Band ',num2str(n)];
        data{n+10,2}=D.Freq_bands{n};
        c=n+10;
end

if isfield(D,'Ridge_recon')
    Ns=size(D.Ridge_recon,1);
    data{c+1,1}='Analysis_type';
    data{c+1,2}=D.Analysis_type;
    data{c+2,1}='Frequency_resolution';
    data{c+2,2}=D.Frequency_resolution;
else         
    Ns=size(D.Filtered_sigs,1); 
end

dstart=c+4;

data{dstart,1}='Time (s)';
for l=1:L;
data{l+dstart,1}=D.Time(l);
end
if ind==1
    if isfield(D,'Filtered_sigs')
    
        for j=1:Ns
            for k=1:Nb
                data{dstart,1+(Ns*(k-1)+j)}=['Sig',num2str(j),' Band',num2str(k)]; 
            for b=1:L
            data{b+dstart,1+(Ns*(k-1)+j)}=D.Filtered_sigs{j,k}(b);
            end
        
            end
    
        end
            
    else  
    for j=1:Ns
        for k=1:Nb
            data{dstart,1+(Ns*(k-1)+j)}=['Sig',num2str(j),' Band',num2str(k)]; 
            for b=1:L
                data{b+dstart,1+(Ns*(k-1)+j)}=D.Ridge_recon{j,k}(b);
            end
        
        end
    
    end
    
    end
elseif ind==2
    
    if isfield(D,'Filtered_sigs')
    
        for j=1:Ns
            for k=1:Nb
                data{dstart,1+(Ns*(k-1)+j)}=['Phase',num2str(j),' Band',num2str(k)]; 
            for b=1:L
            data{b+dstart,1+(Ns*(k-1)+j)}=D.Filtered_phases{j,k}(b);
            end
        
            end
    
        end
            
    else  
    for j=1:Ns
        for k=1:Nb
            data{dstart,1+(Ns*(k-1)+j)}=['Freq',num2str(j),' Band',num2str(k)]; 
            for b=1:L
                data{b+dstart,1+(Ns*(k-1)+j)}=D.Ridge_frequency{j,k}(b);
            end
        
        end
    
    end
    
    end
elseif ind==3
    
    if isfield(D,'Filtered_sigs')
    
        for j=1:Ns
            for k=1:Nb
                data{dstart,1+(Ns*(k-1)+j)}=['Amp',num2str(j),' Band',num2str(k)]; 
            for b=1:L
            data{b+dstart,1+(Ns*(k-1)+j)}=D.Filtered_amplitudes{j,k}(b);
            end
        
            end
    
        end
            
    else  
    for j=1:Ns
        for k=1:Nb
            data{dstart,1+(Ns*(k-1)+j)}=['Amp',num2str(j),' Band',num2str(k)]; 
            for b=1:L
                data{b+dstart,1+(Ns*(k-1)+j)}=D.Ridge_amplitude{j,k}(b);
            end
        
        end
    
    end
    
    end
else
end
