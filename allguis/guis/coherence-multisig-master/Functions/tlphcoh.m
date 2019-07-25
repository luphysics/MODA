%----------------Time-localized wavelet phase coherence--------------------
% TPC = tlphcoh(WT1,WT2,freq,fs,Optional:numcycles)
% calculates time-localized wavelet phase coherence TPC.
%
% Input:
% WT1,WT2 - wavelet transforms of two signals
% freq - frequencies used in wavelet transform
% fs - sampling frequency of a signals from which WT1, WT2 were calculated
% numcycles - number of cycles for calculating TPC (determines adaptive
%             window length, i.e. at 0.1 Hz it will be (1/0.1)*numcycles
%             seconds); default=10.
%
% Author: Dmytro Iatsenko (http://www.physics.lancs.ac.uk/research/nbmphysics/diats)
%--------------------------------------------------------------------------


function TPC = tlphcoh(TFR1,TFR2,freq,fs,varargin)

[NF,L]=size(TFR1);
if nargin>4, wsize=varargin{1}; else wsize=10; end

IPC=exp(1i*angle(TFR1.*conj(TFR2)));
ZPC=IPC; ZPC(isnan(ZPC))=0; cumPC=[zeros(NF,1),cumsum(ZPC,2)];
TPC=zeros(NF,L)*NaN;
for fn=1:NF
    cs=IPC(fn,:); cumcs=cumPC(fn,:);
    tn1=find(~isnan(cs),1,'first'); tn2=find(~isnan(cs),1,'last');
    
    window=round((wsize/freq(fn))*fs); window=window+1-mod(window,2); hw=floor(window/2);
    
    if ~isempty(tn1+tn2) && window<=tn2-tn1
    locpc=abs(cumcs(tn1+window:tn2+1)-cumcs(tn1:tn2-window+1))/window;
    TPC(fn,tn1+hw:tn2-hw)=locpc;
    end
end

end

