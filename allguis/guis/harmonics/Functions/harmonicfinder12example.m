function [outputa,scalefrequency1]=harmonicfinder12example(trans1,detsig,parametri,surrnumber,savedata);

scalefrequency

[m,n]=size(trans1);
res=NaN(m,m);
for a1=1:m
    margin=ceil((sum(isnan(angle(trans1(a1,1:n)))))/2);
    phase1=angle(trans1(a1,1+margin:n-margin)); %slow
    for a2=1:a1
        phase2=angle(trans1(a2,1+margin:n-margin)); %fast
        if a2==a1 && any(phase1-phase2)==1
            error=1
        end
        
        [binner,meanindex]=indexfinder3(phase1,phase2);
        res(a1,a2)=meanindex;
        res(a2,a1)=meanindex;
    end %a2
end %a1

ressurr=NaN(surrnumber,m,m);
for sigb=1:surrnumber
    surrsig=aaft4(detsig');
    transsurr=modbasicwavelet_flow_cmplx4(surrsig, parametri);
    
    sigb
    
    for a1=1:m
        margin=ceil((sum(isnan(angle(trans1(a1,1:n)))))/2);
        phase1=angle(trans1(a1,1+margin:n-margin)); %slow
        
        for a2=1:a1
            phase2=angle(transsurr(a2,1+margin:n-margin)); %fast
            
            [binner,meanindex]=indexfinder3(phase1,phase2);
            ressurr(sigb,a1,a2)=meanindex;
            ressurr(sigb,a2,a1)=meanindex;
            
        end %a2
    end %a1
end %sigb

figure
s = mesh(scalefrequency1,scalefrequency1,res');
s.FaceColor = "flat";

set(gca,'xscale','log')
set(gca,'yscale','log')
title([savedata ' raw harmonics'])
% shading interp; 
view(2)
axis tight
%            [Y,I] = SORT(X,DIM,MODE) also returns an index matrix I.
%     If X is a vector, then Y = X(I).
%     If X is an m-by-n matrix and DIM=1, then
%         for j = 1:n, Y(:,j) = X(I(:,j),j); end
for a1=1:m
    for a2=1:a1
        [ysurr,isurr]=sort([res(a1,a2); ressurr(:,a1,a2)]);
        sig(isurr)=[0:1:surrnumber];
        pos(a1,a2)=sig(1);
        pos(a2,a1)=sig(1);
    end
end
figure
s = mesh(scalefrequency1,scalefrequency1,pos');
s.FaceColor = "flat";

set(gca,'xscale','log')
set(gca,'yscale','log')
title([savedata ' higher than how many aaft surrogates'])
% shading interp; 
view(2)
axis tight

for a1=1:m
    for a2=1:m
        surrmean(a1,a2)=nanmean(ressurr(:,a1,a2));
        surrstd(a1,a2)=nanstd(ressurr(:,a1,a2));
        sig=(res(a1,a2)-surrmean(a1,a2))./surrstd(a1,a2);
        pos(a1,a2)=min(sig, 5); %display maxes out at 5 standard deviations
        %pos(a2,a1)=sig;
    end
end
figure
s = mesh(scalefrequency1,scalefrequency1,pos');
s.FaceColor = "flat";

set(gca,'xscale','log')
set(gca,'yscale','log')
title([savedata ' relative to mean and std of surr distribution'])
% shading interp; 
view(2)
axis tight

outputa=res;
clear transsurr
% save([savedata '.mat'])
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [binner,meanindex]=indexfinder3(data1,data2);

%            [Y,I] = SORT(X,DIM,MODE) also returns an index matrix I.
%     If X is a vector, then Y = X(I).
%     If X is an m-by-n matrix and DIM=1, then
%         for j = 1:n, Y(:,j) = X(I(:,j),j); end
% i is a vector showing where the sorted values came FROM

bins=24; %factor of 2,3,4, and nearly 5
[y1,i1]=sort(data1);
[y2,i2]=sort(data2);
dummy1=[1:1:length(data1)];
dummy2=[1:1:length(data2)];
pslow(i1)=ceil(bins*dummy1/length(data1));
pfast(i2)=ceil(bins*dummy2/length(data2));

% %Initialize RAND to a different state each time.
%            rand('state',sum(100*clock))
%
% %            [Y,I] = SORT(X,DIM,MODE) also returns an index matrix I.
% %     If X is a vector, then Y = X(I).
% %     If X is an m-by-n matrix and DIM=1, then
% %         for j = 1:n, Y(:,j) = X(I(:,j),j); end
%
% amplitudes=randn([1,length(seriesb)]);
%
% [Y1,I1] = sort(amplitudes);
% [Y2,I2] = sort(seriesb);
%
% seriesc(I2)=Y1;
%

binner=zeros(bins,bins);
for n=1:length(pslow)
    binner(pslow(n),pfast(n))=binner(pslow(n),pfast(n))+1;
end
% figure
% mesh(binner)

total=sum(sum(binner));
if total==0
    meanindex=NaN;
else
    i2=0;
    for n=1:bins %this loop finds the entropy i2 of pfast distribution given pslow
        i1(n)=0;
        for m=1:bins %this loop finds the entropy i1 of pfast distribution for each particular bin of pslow
            if sum(binner(n,:))~=0
                pa=(binner(n,m)/sum(binner(n,:)));
            else
                pa=0;
            end
            if pa~=0
                i1(n)=i1(n)-pa*log2(pa);
            end
        end
        pb=sum(binner(n,:))/total;
        i2=i2+(pb*i1(n));
    end
    
    i3=0;
    for n=1:bins %this loop finds the entropy i3 of the pfast distribution without knowing pslow
        if sum(binner(:,n))~=0
            pc=sum(binner(:,n))/sum(sum(binner));
        else
            pc=0;
        end %if
        if pc~=0
            i3=i3-pc*log2(pc);
        end
    end
    
    meanindex=(i3-i2)/i3;
end %if
end

