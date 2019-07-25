function [surr,varargout] = surrogate(signal,N,varargin)
%Calculates Surrogates for a given time series and yields matrix of signals
[m,L]=size(signal);
if L==1, signal=(signal(:))'; [m,L]=size(signal); end %correct if not string

surr=zeros(N,L);

if nargin>2
    method=varargin{1};
else
    method='RP';
end

if strcmp(method,'RP') % Random Permutation surrogates
    for k=1:N
        surr(k,:)=signal(randperm(L));
    end    
    
elseif strcmp(method,'FT') % FT surrogates
    ll=ceil(L/2);
    ftsig=fft(signal,L);
    for k=1:N
        surr(k,1)=ftsig(1); randph=2*pi*randn(1,ll-1);
        surr(k,2:ll)=ftsig(2:ll).*exp(1i*randph);
        surr(k,2+L-ll:L)=conj(fliplr(surr(k,2:ll)));
        
        surr(k,:)=real(ifft(surr(k,:),L));
    end
    
elseif strcmp(method,'AAFT') % AAFT surrogates
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
    
elseif strcmp(method,'IAAFT2') % IAAFT2 (exact spectrum), PS first seed
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
end
end