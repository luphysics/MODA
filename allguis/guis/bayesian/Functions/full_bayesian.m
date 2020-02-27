function [tm,p1,p2,cpl1,cpl2,cf1,cf2,mcf1,mcf2,surr_cpl1,surr_cpl2] = full_bayesian(sig1, sig2, int11, int12, int21, int22, fs, win, pr, ovr, bn, ns, signif)
% Created for use in PyMODA.

int1 = [int11, int12];
int2 = [int21, int22];

[bands1,~] = loop_butter(sig1(:),int1(:),fs);
phi1=angle(hilbert(bands1));

[bands2,~] = loop_butter(sig2(:),int2(:),fs);
phi2=angle(hilbert(bands2));

p1=phi1;
p2=phi2;

[tm, cc, e] = bayes_main(phi1, phi2, win, 1/fs, ovr, pr,0,bn);
for m=1:size(cc,1)
    [cpl1(m),cpl2(m)]=dirc(cc(m,:),bn); % Direction of coupling
    [~,~,q21(:,:,m),q12(:,:,m)]=CFprint(cc(m,:),bn); % Coupling functions
end

cf1 = q21;
cf2 = q12;
mcf1 = squeeze(mean(q21,3));
mcf2 = squeeze(mean(q12,3));

surr1 = surrcalc(phi1',ns,'CPP',0,fs);
surr2 = surrcalc(phi2',ns,'CPP',0,fs);

for n=1:ns
    [~,cc_surr{n}]=bayes_main(surr1(n,:),surr2(n,:),win,1/fs,ovr,pr,1,bn);
    for idx=1:size(cc_surr{n},1)
        [scpl1(n,idx),scpl2(n,idx)]=dirc(cc_surr{n}(idx,:),bn);
    end
end

alph=signif;
alph=(1-(alph/100));
if floor((ns+1)*alph)==0
    surr_cpl1 = max(scpl1);
    surr_cpl2 = max(scpl2);
else
    K=floor((ns+1)*alph);
    s1=sort(scpl1,'descend');
    s2=sort(scpl2,'descend');
    surr_cpl1 = s1(K,:);
    surr_cpl2 = s2(K,:);
    
end

end

