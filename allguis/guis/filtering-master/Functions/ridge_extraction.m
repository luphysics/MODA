function [transform, freq, iamp, iphi, ifreq, recon] = ridge_extraction(calc_type, signal, fs, interval1, interval2, cut_edges, preprocess, w_type)


if calc_type == 1
    [WT,freqarr,wopt]=wt(signal,fs,'fmin',interval1,'fmax',interval2,'CutEdges','off',...
        'Preprocess',preprocess,'Wavelet',w_type);
else
    [WT,freqarr,wopt]=wft(signal,fs,'fmin',interval1,'fmax',interval2,'CutEdges','off',...
        'Preprocess',preprocess,'Window',w_type);
end

tfsupp = ecurve(WT,freqarr,wopt);

[iamp,iphi,ifreq] = rectfr(tfsupp,WT,freqarr,wopt,'direct');

recon = iamp.*cos(iphi);
iphi = mod(iphi,2*pi);

transform = WT;
freq = freqarr;


