
function [tm,cc,e]=bayes_main(ph1,ph2,win,h,ovr,pr,s,bn)
%the main function for the inference, propagation and other evaluations

%---inputs---
%ph1,ph2 - phase time-series vectors
%win - window in seconds
%h - sampling step e.g. h=0.01
%ovr - overlaping of windows; ovr=1 is no overlap; ovr=0.75 will overlap
%      the last 1/4 of each window with the next window
%pr - propagation constant
%s - print progress status if s=1 
%bn - order of Fourier base function

%---outputs---
%tm - time vector for plotting 
%cc - inferred mean parameters
%e  - inferred noise

%example for default input parameters and call of 
%the function >> [tm,cc,e]=bayes_main(ph1,ph2,40,0.01,1,0.2,1,2);
%%

cc=[]; 
win=win/h;
w=ovr*win;
ps=ph2-ph1;
pw=win*h*pr;

M=2+2*((2*bn+1)^2-1);
L=2;

Cpr=zeros(M/L,L);XIpr=zeros(M);


%unwrap the phases if they are not 
if (max(ph1)<(2*pi+0.1))
    ph1=unwrap(ph1);
    ph2=unwrap(ph2);
end

%set the right dimensions for the vectors
[m,n]=size(ph1);
if m<n
    ph1=ph1';
    ph2=ph2';
end


%% do the main calculations for each window
for i=0:floor((length(ps)-win)/w)
    
    phi1=ph1(i*w+1:i*w+win); phi2=ph2(i*w+1:i*w+win);
    
    %-----bayesian inference for one window------
    [Cpt,XIpt,E]=bayesPhs(Cpr,XIpr,h,500,0.00001,phi1',phi2',bn);
   
   
     %the propagation for the next window
     [XIpr,Cpr] = Propagation_function_XIpt(Cpt,XIpt,pw);
    
  
 
 e(i+1,:,:)=E; 
 cc(i+1,:)=Cpt(:);
   
 %display progress 
 if s 
   display(['processed so far: t= ' num2str((i+1)*w*h) 's /' num2str(length(ph1)*h) 's ;']);
   %Cpt
   %E
 end
 
end

%time vector for plotting
tm = (win/2:w:length(ph1)-win/2)*h;
%%

 

function [ XIpr,Cpr ] = Propagation_function_XIpt(Cpt,XIpt,p)
% Propagation function with covariance
% find the new prior for the next block

% The average is not supposed to change
Cpr=Cpt;   


% Prepare the diffusion matrix
% Set a value for a particular parameter
Inv_Diffusion=zeros(length(XIpt));
invXIpt=inv(XIpt);
for i=1:length(Cpt(:)) 
    Inv_Diffusion(i,i) = p*p* invXIpt(i,i);
end

% The gaussian of the posterior is convoluted with another 
% gaussian which express the diffusion of the parameter.
XIpr=inv(( invXIpt + Inv_Diffusion ));
%%

 function [ XIpr,Cpr ] = Propagation_function_Cpt(Cpt,XIpt,p)
% Propagation function with parameters
% find the new prior for the next block

% The average is not supposed to change
Cpr=Cpt;   


% Prepare the diffusion matrix
% Set a value for a particular parameter
Inv_Diffusion=zeros(length(XIpt));
invXIpt=inv(XIpt);
for i=1:length(Cpt(:)) 
    Inv_Diffusion(i,i) = p*p* (Cpt(i));
end

% The gaussian of the posterior is convoluted with another 
% gaussian which express the diffusion of the parameter.
XIpr=inv(( invXIpt + Inv_Diffusion ));
%%