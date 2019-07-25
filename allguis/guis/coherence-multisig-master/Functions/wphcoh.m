% [phcoh,Optional:phdiff] = wphcoh(WT1,WT2)
% returns time-averaged wavelet phase coherence between two signals;
% WT1 and WT2 are wavelet transforms of these signals.
%
% Author: Dmytro Iatsenko (http://www.physics.lancs.ac.uk/research/nbmphysics/diats)
%--------------------------------------------------------------------------

function [phcoh,varargout] = wphcoh(WT1,WT2)

FN=min([size(WT1,1),size(WT2,1)]);
WT1=WT1(1:FN,:); WT2=WT2(1:FN,:);
phi1=angle(WT1); phi2=angle(WT2);
phexp=exp(1i*(phi1-phi2));

phcoh=zeros(1,FN)*NaN; phdiff=zeros(1,FN)*NaN;
for fn=1:FN
    cphexp=phexp(fn,:); cphexp=cphexp(~isnan(cphexp));
    NL=length(find(WT1(fn,:)==0 & WT2(fn,:)==0));
    CL=length(cphexp);
    if CL>0
        phph=mean(cphexp)-NL/CL;
        phcoh(fn)=abs(phph);
        phdiff(fn)=angle(phph);
    end
end

if nargout>1, varargout{1}=phdiff; end

end

