
scale_min = parametri.valcki.scale_min;
scale_max = parametri.valcki.scale_max;
sigma = parametri.valcki.sigma;

clear scalefrequency1

m_max = floor(log(scale_max/scale_min)/log(sigma));
m = 0:1:m_max+1;

for z = 1 : length(m)
    
    scalefrequency1(z)=1/(scale_min*(1.05^(z-1)));
    
end

clear scale_min scale_max sigma m_max m z