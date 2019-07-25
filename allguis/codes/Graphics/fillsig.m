% Fill between curves based on threshold

function val=fillsig(xvec,sig1,sig2,thresh1,thresh2,direction,col)
% sig1, sig2 - signals to shade between
% thresh1, thresh2 - signals to calculate significance
% xvec=freq;
% sig1=pLANH_PNH;
% sig2=x;
% direction='less'; % Can be 'greater'
% col=[.5 .5 .5];
% 
% figure
% semilogx(xvec,sig1)
% hold on
% semilogx(xvec,sig2,'r--');

if strcmp(direction,'less')
    n=find(thresh1<thresh2);
elseif strcmp(direction,'greater')
    n=find(thresh1>thresh2);
end

if isempty(n)
    val=0;
else
    val=1;
end

d=diff(n);
n2=find(d>1);

if isfinite(n2)

if length(n2)==1
    ar{1}=n(1:n2);
    ar{2}=n(n2+1:end);
else
    ar{1}=n(1:n2(1));
    ar{length(n2)+1}=n(n2(end)+1:end);
    for k=1:length(n2)-1
        ar{k+1}=n(n2(k)+1:n2(k+1));
    end
        
end
else
    ar{1}=n;
end
s=size(ar);
for j=1:s(2)
    f=xvec(ar{j});
    X=[f',fliplr(f')];
    Y=[sig1(ar{j}),fliplr(sig2(ar{j}))];
    fill(X,Y,col,'EdgeColor',col)
end

