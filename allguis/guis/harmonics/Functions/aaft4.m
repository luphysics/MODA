function [output,finalorder]=aaft4(seriesa);

[idim jdim] = size(seriesa);
if idim > 1
    seriesa = seriesa';
end
%Initialize RAND to a different state each time.
rand('state',sum(100*clock))

%            [Y,I] = SORT(X,DIM,MODE) also returns an index matrix I.
%     If X is a vector, then Y = X(I).
%     If X is an m-by-n matrix and DIM=1, then
%         for j = 1:n, Y(:,j) = X(I(:,j),j); end
% i describes where each value of y came FROM in x

amplitudes=randn([1,length(seriesa)]);
[Y1,I1] = sort(amplitudes);
[Y2,I2] = sort(seriesa);
seriesc(I2)=Y1;


fftseries=fft(seriesc);


if mod(length(seriesa),2)==1
    
    numbers=rand([1,((length(fftseries)-1)/2)]);
    phases=2*pi*numbers;
    phasors=exp(i*phases);
    
    % test1=length(fftseries)
    % test2=length([1, phasors, 1, fliplr(conj(phasors))])
    
    rotfftseries=fftseries.*[1, phasors, fliplr(conj(phasors))];
    
    ifftrot=ifft(rotfftseries);
    
    [Y3,I3] = sort(ifftrot);
    
    output(I3)=Y2;
else
    
    numbers=rand([1,(length(fftseries)/2-1)]);
    phases=2*pi*numbers;
    phasors=exp(i*phases);
    
    % test1=length(fftseries)
    % test2=length([1, phasors, 1, fliplr(conj(phasors))])
    
    rotfftseries=fftseries.*[1, phasors, 1, fliplr(conj(phasors))];
    
    ifftrot=ifft(rotfftseries);
    
    [Y3,I3] = sort(ifftrot);
    
    output(I3)=Y2;
end %if
origorder=[1:1:length(seriesa)];
sortedorder=origorder(I2);
finalorder(I3)=sortedorder;
end