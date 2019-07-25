%funcion generating surrogates: read further for different methods, 
%most common is random permutation 'RP' or phase shuffle 'FT'.
%N= munber of surrogates to generate


function [surr,varargout] = surrogate(signal,N,varargin)
[m,L]=size(signal);

surr=zeros(N,L);

if nargin>2
    method=varargin{1};
else
    method='RP';
end

if strcmp(method,'RP') % Random Permutation surrogates
    %New data are created simply by random permutations of the original series.
%     The permutations guarantee the same amplitude distribution than the original series, 
%     but destroy any linear correlation. This method is associated to the null hypothesis of the 
%     data being uncorrelated noise (possibly Gaussian and measured by a static nonlinear function
    for k=1:N
        surr(k,:)=signal(randperm(L));
    end    
    
elseif strcmp(method,'FT') % FT surrogates
%     Random Phases; also known as FT, for Fourier Transform):
%         In order to preserve the linear correlation (the periodogram) of the series, 
%         surrogate data are created by the inverse Fourier Transform of the modules of Fourier Transform 
%         of the original data with new (uniformly random) phases. If the surrogates must be real, 
%         the Fourier phases must be antisymmetric with respect to the central value of data.
    ll=ceil(L/2);
    ftsig=fft(signal,L);
    for k=1:N
        surr(k,1)=ftsig(1); randph=2*pi*randgen(1,ll-1);
        surr(k,2:ll)=ftsig(2:ll).*exp(1i*randph);
        surr(k,2+L-ll:L)=conj(fliplr(surr(k,2:ll)));
        
        surr(k,:)=real(ifft(surr(k,:),L));
        endnono
    end
elseif strcmp(method,'AAFT') % AAFT surrogates Amplitude Adjusted Fourier Transform):
%     This method has approximately the advantages of the two previous ones: it tries to 
%     preserve both the linear structure and the amplitude distribution.

    ll=ceil(L/2);
    sigma=std(signal);
    [sortsig,sortind]=sort(signal); rankind(sortind)=1:L;
    rescrank=zeros(1,L);
    for k=1:N
        rgs=sort(randgen(1,L,0,sigma));
        rescsig=rgs(rankind);
        
        ftresc=fft(rescsig,L);
        surr(k,1)=ftresc(1); randph=2*pi*randgen(1,ll-1);
        %surr(k,2:ll)=abs(ftresc(2:ll)).*exp(1i*randph); % should we randomize all phases, since there might be some bias?
        surr(k,2:ll)=ftresc(2:ll).*exp(1i*randph);
        surr(k,2+L-ll:L)=conj(fliplr(surr(k,2:ll)));
        
        surr(k,:)=real(ifft(surr(k,:),L));
        
        [~,rescind]=sort(surr(k,:)); rescrank(rescind)=1:L;
        surr(k,:)=sortsig(rescrank);
    end
    
elseif strcmp(method,'AAFTFT') % AAFT with FT distribution instead of noise
    ll=ceil(L/2);
    ftsig=fft(signal,L);
    [sortsig,sortind]=sort(signal); rankind(sortind)=1:L;
    rescrank=zeros(1,L);
    for k=1:N
        surr(k,1)=ftsig(1); randph=2*pi*randgen(1,ll-1);
        surr(k,2:ll)=ftsig(2:ll).*exp(1i*randph);
        surr(k,2+L-ll:L)=conj(fliplr(surr(k,2:ll)));
        surr(k,:)=real(ifft(surr(k,:),L));
        
        rgs=sort(surr(k,:));
        rescsig=rgs(rankind);
        
        ftresc=fft(rescsig,L);
        surr(k,1)=ftresc(1); randph=2*pi*randgen(1,ll-1);
        surr(k,2:ll)=ftresc(2:ll).*exp(1i*randph);
        surr(k,2+L-ll:L)=conj(fliplr(surr(k,2:ll)));
        
        surr(k,:)=real(ifft(surr(k,:),L));
        
        [~,rescind]=sort(surr(k,:)); rescrank(rescind)=1:L;
        surr(k,:)=sortsig(rescrank);
    end
    
elseif strcmp(method,'IAAFT') % IAAFT2 (exact spectrum), PS first seed
%      Iterative Amplitude Adjusted Fourier Transform): This algorithm is an iterative version of AAFT.
%      The steps are repeated until the autocorrelation function is sufficiently similar to the original, 
%      or until there is no change in the amplitudes.
    maxitn=1000; % maximum number of iterations
    [sortsig,sortind]=sort(signal); rankind(sortind)=1:L;
    ftsig=fft(signal,L);
    ovitn=0;
    for  k=1:N
        surr(k,:)=signal(randperm(L));       
        
        itn=1; iterrank=rankind; olditrank=zeros(1,L);
        while (max(abs(olditrank-iterrank))~=0 & itn<maxitn)
            olditrank=iterrank;
            iterf=real(ifft(abs(ftsig).*exp(1i*angle(fft(surr(k,:),L))))); % replace Fourier amplitudes (real() since makes mistakes of order \epsilon)
            [~,iterind]=sort(iterf); iterrank(iterind)=1:L;
            surr(k,:)=sortsig(iterrank);
            itn=itn+1;
        end
        ovitn=ovitn+itn;
        surr(k,:)=iterf;
    end
    ovitn=ovitn/N;
    
elseif strcmp(method,'IAAFT1') % IAAFT1 (exact distribution), PS first seed
    maxitn=1000; % maximum number of iterations
    [sortsig,sortind]=sort(signal); rankind(sortind)=1:L;
    ftsig=fft(signal,L);
    ovitn=0;
    for  k=1:N
        surr(k,:)=signal(randperm(L));        
        
        itn=1; iterrank=rankind; olditrank=zeros(1,L);
        while (max(abs(olditrank-iterrank))~=0 & itn<maxitn)
            olditrank=iterrank;
            iterf=real(ifft(abs(ftsig).*exp(1i*angle(fft(surr(k,:)))))); % replace Fourier amplitudes
            [~,iterind]=sort(iterf); iterrank(iterind)=1:L;
            surr(k,:)=sortsig(iterrank);
            itn=itn+1;
        end
        ovitn=ovitn+itn;
    end
    ovitn=ovitn/N; % overall number of iterations
    
elseif strcmp(method,'TS') % TS: twin surrogates
%     %this technique generates surrogates which correspond to an inde-pendent copy 
%     of the underlying system, i. e. they induce a trajectory of the underlying 
%     systemstarting at di?erent initial conditions. We show that these surrogates are
%     well suited to test forcomplex synchronisation
    if nargin>3
        dL=varargin{2};
    else
        dL=L;
    end
    alpha=0.1;
    
    Rij=zeros(L,L);
    for k=2:L
        Rij(k,1:k-1)=max(abs(signal(:,1:k-1)-signal(:,k)*ones(1,k-1)));
    end
    Rij=Rij+Rij';
    [~,pl]=min(Rij(1:round(L/2),L));
    Sij=sort(Rij(:)); delta=Sij(round(alpha*L^2)); clear Sij;
    Rij(Rij<delta)=-1; Rij(Rij>delta)=0; Rij=abs(Rij);
    
    ind=cell(L,1); eln=zeros(L,1); twind=1:L;
    remp=1; % remaining points
    while ~isempty(remp)
        twn=remp(1);
        ind{twn}=remp(max(abs(Rij(:,remp)-Rij(:,twn)*ones(1,numel(remp))))==0);
        ind(ind{twn})=ind(twn);
        eln(ind{twn})=length(ind{twn});
        twind(ind{twn})=0;
        remp=twind(twind>0);
    end
    clear Rij twind;
    
    for sn=1:N
        kn=randi(L,1)-1;
        for j=1:dL
            kn=kn+1;
            surr(sn,j)=signal(1,kn);
            kn=ind{kn}(randi(eln(kn),1));
            if kn==L
                kn=pl;
            end
        end
    end
    
elseif strcmp(method,'PPS') % PPS surrogates
%     can distinguish between a noisy periodic
% orbit and deterministic non-periodic inter-cycle dynamics. Possible origins of deterministic nonperiodic
% inter-cycle dynamics include: non-periodic linear or nonlinear dynamics, or chaos. This
% new algorithm is based on mimicking the large-scale dynamics with a local model, but obliterating
% the fine scale features with dynamic noise. 
    if nargin>3
        dL=varargin{2};
    else
        dL=L;
    end
    ssig=zeros(1,L); mind=zeros(1,L);
    for k=1:L
        matr=max(abs(signal(:,:)-signal(:,k)*ones(1,L)));
        [ssig(k),mind(k)]=min(matr(matr>0));
    end
    [~,pl]=min(matr(1:round(L/2))); rho=0.7*mean(ssig); clear mind ssig;
    
    if m==1
    for sn=1:N
        kn=randi(L,1)-1;
        for j=1:dL
            kn=kn+1;
            surr(sn,j)=signal(1,kn);
            sigdist=abs(signal-(signal(kn)+randgen(1,1,0,rho)));
            [~,kn]=min(sigdist);
            if kn==L
                kn=pl;
            end
        end
    end
    else
    %bt=0; flag=0;
    for sn=1:N
        kn=randi(L,1)-1;
        for j=1:dL
            kn=kn+1;
            %knold=kn;
            surr(sn,j)=signal(1,kn);
            sigdist=max(abs(signal(:,:)-(signal(:,kn)+randgen(m,1,0,rho))*ones(1,L)));
            [~,kn]=min(sigdist);
            if kn==L
                kn=pl;
            end
            %{
            if knold==kn
                flag=flag+1;
            else
                if flag==1
                    bt=bt+1;
                end
                flag=0;
            end
            %}
        end
    end
    
    end
    
elseif strcmp(method,'CPP') % cycle phase permutation surrogates
    
    signal=mod(signal,2*pi);
    
    dcpoints=find(signal(2:end)-signal(1:end-1)<-pi);
    NC=length(dcpoints)-1;
    if NC>0
    cycles=cell(NC,1);
    for k=1:NC
        cycles{k}=signal(dcpoints(k)+1:dcpoints(k+1));
    end
    stcycle=signal(1:dcpoints(1));
    endcycle=signal(dcpoints(k+1)+1:end);
    
    for sn=1:N
        surr(sn,:)=unwrap(horzcat(stcycle,cycles{randperm(NC)},endcycle));
    end
    
    else
    for sn=1:N
        surr(sn,:)=unwrap(signal);
    end
    end
    
elseif strcmp(method,'MCPP') % cycle phase permutation surrogates
    
    signal=mod(signal,2*pi);
    
    dcpoints=find(signal(2:end)-signal(1:end-1)<-pi);
    fcpoints=find(signal(2:end)-signal(1:end-1)>pi);
    NC=length(dcpoints)-1;
    if NC>1
    cycles=cell(NC+1,1);
    cycles{1}=signal(1:dcpoints(1));
    cn=1;
    for k=1:NC
        if isempty(fcpoints(fcpoints>dcpoints(k) & fcpoints<dcpoints(k+1)))
            cn=cn+1;
            cycles{cn}=signal(dcpoints(k)+1:dcpoints(k+1));
        else
            cycles{cn}=horzcat(cycles{cn},signal(dcpoints(k)+1:dcpoints(k+1)));
        end
    end
    stcycle=cycles{1};
    endcycle=signal(dcpoints(end)+1:end);
    cycles=cycles(2:cn);
    for sn=1:N
        surr(sn,:)=unwrap(horzcat(stcycle,cycles{randperm(cn-1)},endcycle));
    end
    
    else
    for sn=1:N
        surr(sn,:)=unwrap(signal);
    end
    end
    
elseif strcmp(method,'tshift') % time shift permutation surrogates
    
    for sn=1:N
        if nargin>3
           startp=varargin{2};
        else
           startp=randi(L-1,1);
        end
        surr(sn,:)=horzcat(signal(1+startp:L),signal(1:startp));
        if nargout>1
            varargout{1}(sn)=startp;
        end
    end
    
elseif strcmp(method,'tshift2') % time shift permutation surrogates
    
    cutp=ceil(L/2);
    
    surr=zeros(N,cutp);
    for sn=1:N
        if nargin>3
            startp=varargin{2};
        else
            startp=randi(L-cutp,1);
        end
        surr(sn,:)=signal(1+startp:startp+cutp);
    end
    
elseif strcmp(method,'CAAFT') % CAAFT (or STAP with m=1)
    ll=ceil(L/2);
    ms=mean(signal); sigma=std(signal);
    [sortsig,sortind]=sort(signal); rankind(sortind)=1:L;
    rescrank=zeros(1,L);
        
    tmax=floor(L/4); p=tmax; K=N; % parameters of CAAFT
    y=zeros(1,L); yft=zeros(1,L); z=zeros(1,L);
    rx=zeros(1,tmax+1); ry=zeros(1,tmax+1); rz=zeros(1,tmax+1); % autocorrelations
    ar=zeros(1,L); arrank=zeros(1,L); arcf=zeros(K,p+1); arstd=zeros(K,1); rw=zeros(K,tmax+1);
    
    for tn=0:tmax % let us take without normalization
        rx(tn+1)=mean((signal(1:L-tn)-ms).*(signal(1+tn:L)-ms));
    end
    
    for k=1:K
        rgs=sort(randgen(1,L,0,sigma));
        y=rgs(rankind);
        
        ftresc=fft(y,L);
        yft(1)=ftresc(1); randph=2*pi*randgen(1,ll-1);
        yft(2:ll)=ftresc(2:ll).*exp(1i*randph);
        yft(2+L-ll:L)=conj(fliplr(yft(2:ll)));
        
        yft=ifft(yft,L);
        
        [~,rescind]=sort(yft); rescrank(rescind)=1:L;
        z=sortsig(rescrank);
        
        my=mean(y); mz=mean(z);
        for tn=0:tmax % let us take without normalization
            ry(tn+1)=mean((y(1:L-tn)-my).*(y(1+tn:L)-my));
            rz(tn+1)=mean((z(1:L-tn)-mz).*(z(1+tn:L)-mz));
        end
        
        cf=polyfit(rz,ry,1);
        ru=cf(2)+cf(1)*rx;
        
        [arcf(k,:),arstd(k)]=levinson(ru,p); arstd(k)=sqrt(arstd(k)); arcf(k,:)=-arcf(k,:);
        for fn=1:p
            for kn=1:fn-1
                ar(fn)=ar(fn)+arcf(k,kn+1)*ar(fn-kn);
            end
            ar(fn)=ar(fn)+randgen(1,1,0,arstd(k));
        end
        for arn=1+p:L
            ar(arn)=ar(arn-p:arn-1)*flipud(arcf(k,2:end)')+randgen(1,1,0,arstd(k));
        end
        
        [~,indar]=sort(ar); arrank(indar)=1:L;
        w=sortsig(arrank); mw=mean(w);
        
        for tn=0:tmax
            rw(k,tn+1)=mean((w(1:L-tn)-mw).*(w(1+tn:L)-mw));
        end
    end
    
    [~,indm]=min(std(rw-ones(K,1)*rx,0,2));
    arcf=arcf(indm,:); arstd=arstd(indm); rw=rw(indm,:);
    for k=1:N
        for fn=1:p
            for kn=1:fn-1
                ar(fn)=ar(fn)+arcf(kn+1)*ar(fn-kn);
            end
            ar(fn)=ar(fn)+randgen(1,1,0,arstd);
        end
        for arn=1+p:L
            ar(arn)=ar(arn-p:arn-1)*flipud(arcf(2:end)')+randgen(1,1,0,arstd);
        end
        [~,indar]=sort(ar); arrank(indar)=1:L;
        surr(k,:)=signal(arrank);
    end
    
end
    
end