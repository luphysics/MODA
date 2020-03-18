close all
clear all

time=100
freq=(1/7)

sampl_freq=50;

times=[1/sampl_freq:1/sampl_freq:time];
% t_series=0.1*randn(1,length(times))+sin(2*pi*freq*times).^3;
load("t_series.mat");

figure; plot(times,t_series)
xlabel('times');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale_min=0.5;
scale_max=40;
sigma=1.05;
time_res=0.1;

parametri.valcki.sig_sampl_freq=sampl_freq;
parametri.valcki.scale_min=scale_min;
parametri.valcki.scale_max= scale_max;
parametri.valcki.sigma=sigma;
parametri.valcki.obs_time_res=time_res;

scalefrequency
[output] = modbasicwavelet_flow_cmplx4(t_series, parametri);

figure; surf([0:time_res:time],1./scalefrequency1,abs(output)); shading interp
xlabel('times'); ylabel('1/freq')

[outputa,scalefrequency1]=harmonicfinder12example(output,t_series,parametri,4,'test250220');
%set(gca,'xtick',[1/7 2/7 3/7],'xticklabel',[1/7 2/7 3/7],'FontSize',6);
%set(gca,'ytick',[1/7 2/7 3/7],'yticklabel',[1/7 2/7 3/7],'FontSize',6);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
