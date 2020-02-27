%%% the output is a complex value of the transformation %%%

function [output] = modbasicwavelet_flow_cmplx4(t_series, parametri)

% parametri.valcki.sig_sampl_freq = 10;
% parametri.valcki.sig_t_start = 0;
% parametri.valcki.sig_t_fin = length(det1)/10;
% parametri.valcki.scale_min = 0.4;
% parametri.valcki.scale_max = 400;
% parametri.valcki.sigma = 1.05;
% parametri.valcki.obs_time1 = 0;
% parametri.valcki.obs_time2 = length(det1)/10;
% %there is no sample at t=0
% parametri.valcki.obs_time_res = 1;

sampl_freq = parametri.valcki.sig_sampl_freq;

scale_min = parametri.valcki.scale_min;
scale_max = parametri.valcki.scale_max;
sigma = parametri.valcki.sigma;

time_res = parametri.valcki.obs_time_res;

t_start = 0;
t_fin = length(t_series)/sampl_freq;


m_max = floor(log(scale_max/scale_min)/log(sigma));
m = 0:1:m_max+1;

REZ = NaN(length(m),  ((floor(t_fin - t_start))/time_res) + 1);

flo=floor((t_fin - t_start));
stevec = 0;

for z = 1 : length(m)
    
    s = scale_min*sigma.^m(z);
    
    %begin calculating wavelet
    
    f0 = 1;
    
    tval_k=0;
    tval_k2=0;
    
    for ttt = 0 : 1/sampl_freq : 5000
        zzz = s^(-0.5)*exp(-ttt*ttt/(2*s*s));
        if zzz < 0.1*(s^(-0.5)) && tval_k2==0
            tval_k2 = ttt;
            check1(z)=ttt;
            
        end
        if zzz < 0.001*(s^(-0.5))
            tval_k = ttt  ;  % cas trajanja valcka
            check2(z)=ttt;
            break
        end
    end
    
    st_kor = (tval_k)*sampl_freq;
    margin = (tval_k2)*sampl_freq;
    %round up st_kor for accuracy
    if mod(st_kor,sampl_freq) ~= 0
        st_kor = sampl_freq*ceil(st_kor/sampl_freq);
    end
    %round up margin for safety
    if mod(margin,sampl_freq) ~= 0
        margin = sampl_freq*ceil(margin/sampl_freq);
    end
    %margin determines how close large wavelets can come to the edges of the
    %timeseries
    if margin > 0.5*sampl_freq*flo
        break
    end
    
    u = -st_kor/sampl_freq : 1/sampl_freq : st_kor/sampl_freq;
    
    wavelet = s^(-0.5)*exp(-i*2*pi*f0.*u/s) .* exp(-u.*u/(2*s*s)) ;
    
    X = fft([wavelet zeros(1,length(t_series)-1)]);
    Y = fft([t_series zeros(1,length(wavelet)-1)]);
    con = ifft(X.*Y);
    %con = conv(t_series,wavelet);
    
    rez = con(st_kor + margin + 1 : sampl_freq*time_res : (st_kor - margin + (flo*sampl_freq) +1));
    
    % margin
    % length(rez)
    %rez2= rez(margin + 1 : sampl_freq : (length(rez)-margin)+1);
    stevec = stevec + 1;
    if margin/(sampl_freq*time_res) > 0
        trans = [NaN(1,(margin/(sampl_freq*time_res))) rez NaN(1,(margin/(sampl_freq*time_res)))];
    else
        trans = rez;
    end
    lengthtrans = length(trans);
    lengthrez = length(REZ(stevec,:));
    if lengthtrans==lengthrez
        REZ(stevec,:) = trans;
    end
    
end
% % check1
% % check2
output = [REZ];
end