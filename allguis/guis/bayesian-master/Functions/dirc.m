function [cpl1,cpl2,drc]=dirc(c,bn)   
%calculates the net couplings as norms from the relevant inferred parameters 

%---inputs---
%c - vector of inferred parameters
%bn - order of Fourier base function

%---outputs---
%cpl1 - coupling from second to first oscillator
%cpl2 - coupling from first to second oscillator
%drc  - direction of coupling drc~[-1,1]

%Note that the input is vector of parameters for one time window - for all
%time windows use 'dirc.m' in loop; see e.g. 'example2_CplPncrPrm.m'

%% ------------------------------------------------------------------------                    
                    q1=[]; q2=[];
                    iq1=1; iq2=1; 
                    br=2;
                    K=length(c)/2;
                   
                    for ii=1:bn 
                        q1(iq1)=c(br);     iq1=iq1+1; q1(iq1)=c(br+1);  iq1=iq1+1;
                        q2(iq2)=c(K+br); iq2=iq2+1; q2(iq2)=c(K+br+1);  iq2=iq2+1;
                        br=br+2;
                    end
                    
                    for ii=1:bn 
                        q1(iq1)=c(br);     iq1=iq1+1; q1(iq1)=c(br+1);  iq1=iq1+1;
                        q2(iq2)=c(K+br); iq2=iq2+1; q2(iq2)=c(K+br+1);  iq2=iq2+1;
                        br=br+2;                        
                    end
                    
                    for ii=1:bn
                        for jj=1:bn        
                                    q1(iq1)=c(br);     iq1=iq1+1; q1(iq1)=c(br+1);  iq1=iq1+1;
                                    q2(iq2)=c(K+br); iq2=iq2+1; q2(iq2)=c(K+br+1);  iq2=iq2+1;
                                    br=br+2;
                                   
                                     q1(iq1)=c(br);     iq1=iq1+1; q1(iq1)=c(br+1);  iq1=iq1+1;
                                     q2(iq2)=c(K+br); iq2=iq2+1; q2(iq2)=c(K+br+1);  iq2=iq2+1;
                                     br=br+2;
                        end
                    end                    

            
            cpl1=norm(q1);
            cpl2=norm(q2);
            drc=(cpl2-cpl1)/(cpl1+cpl2);
