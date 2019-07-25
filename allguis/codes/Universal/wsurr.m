% WIAAFT surrogate
clear surrLev
sig=y;
%wvlt='la16';
fs=1;

w=modwt(sig);
N2=size(w,1);
N=length(sig);

for j=1:N2
    tmp=surrcalc(sig, 1, 'IAAFT2', 0, fs);
    surrLev(j,:)=matchRotation(w(j,:),tmp);
    
end

surr=imodwt(surrLev);

maxit = 200; % maximum number of iterations
    [sorted, sortInd] = sort(sig); 
    itrank(sortInd) = linspace(1, N, N); 
    ftsig = fft(sig);
    surrtmp = surr;       
    iter = 1; 
    oldrank = zeros(1, N);

    while (max(abs(oldrank - itrank)) ~= 0 && iter < maxit) % equal spectrum, similar amplitude distribution
        oldrank = itrank;
        % replace Fourier amplitudes (real() since makes mistakes of order \epsilon)
        itFt = real(ifft(abs(ftsig).* exp(1i * angle(fft(surrtmp))))); 
        [~, itind] = sort(itFt); % find the ordering of the new signal 
        itrank(itind)= linspace(1, N, N); 
        surrtmp = sorted(itrank); 
        iter = iter + 1;
    end
    surr = itFt;
    
    figure
    plot(time,sig)
    hold on
    plot(time,surr)
    legend('Sig','Surr')