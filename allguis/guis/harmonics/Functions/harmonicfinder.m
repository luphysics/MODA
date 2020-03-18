function [outputArg1,outputArg2] = harmonicfinder(signal, fs, scale_min, scale_max, sigma, time_res, surr_count)

    scalefreq = scale_frequency(scale_min, scale_max, sigma)
    [output] = modbasicwavelet_flow_cmplx4(signal, );

end

