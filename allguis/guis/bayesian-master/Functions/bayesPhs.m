function [Cpt, XIpt, E]=bayesPhs(Cpr,XIpr,h,max_loops,eps,phi1,phi2,bn)
%infers the parameters and the noise within one block of data

%---inputs---
%Cpr - prior vector of parameters
%XIpr - prior concentration matrix (i.e. inv(covariance) matrix)
%h - sampling rate (e.g. h=0.01)
%max_loops - number of loops to stop for the error for convergence
%eps - error/tolerance for convergence
%phi1,phi2 - input block of data
%bn - order of Fourier base function

%---outputs---
%Cpt - posterior vector of parameters
%XIpt - poaterior concentration matrix (i.e. inv(covariance) matrix)
%E - noise matrix

%% ------------

%phi1S & phi2S - the midpoint time series
phi1S=(phi1(2:end)+phi1(1:end-1))/2; phi2S=(phi2(2:end)+phi2(1:end-1))/2;

%phi1T & phi2T - the derivative time series
phi1T=(phi1(2:end)-phi1(1:end-1))/h; phi2T=(phi2(2:end)-phi2(1:end-1))/h;
phiT=[phi1T;phi2T]; 
 
L=2; %number of dimensions 
M=2+2*((2*bn+1)^2-1); %number of base functions
K=M/L; %number of dimensions for the C vector
lastwarn(''); %reset the last warning message

%the following evaluations fill temporary variables to be used bellow in the main calculations 
p=calculateP(phi1S,phi2S,K,bn);
v1=calculateV(phi1S,phi2S,K,bn,1);
v2=calculateV(phi1S,phi2S,K,bn,2);


%loop to converge iteratively the result until desired precision is reached
C_old=Cpr; 
Cpt=Cpr; %initialize Cpt
for loop=1:max_loops
   
    [E]=calculateE(Cpt',phiT,L,h,p);    
    [Cpt,XIpt]=calculateC(E,p,v1,v2,Cpr,XIpr,M,L,phiT,h);
        
    if(sum((C_old-Cpt).*(C_old-Cpt)./(Cpt.*Cpt))<eps)      
		return; 
    end
     C_old=Cpt;
end
display 'Warning: Desired precision not reached after max_loops.'  
%%

%----------additional functions: -------

%evaluates the base functions 
function [p]=calculateP(phi1,phi2,K,bn)   
    p(K,length(phi1))=0;
    
    p(1,:)=1;    
    br=1;

    for i=1:bn
                p(br+1,:)    = sin(i*phi1);
                p(br+2,:)    = cos(i*phi1);
                br=br+2;
    end
      
    for i=1:bn
                p(br+1,:)    = sin(i*phi2);
                p(br+2,:)    = cos(i*phi2);
                br=br+2;                
    end

    for i=1:bn
        for j=1:bn       
                p(br+1,:)    = sin(i*phi1+j*phi2);                
                p(br+2,:)    = cos(i*phi1+j*phi2);
                br=br+2;                 
                p(br+1,:)    = sin(i*phi1-j*phi2);                
                p(br+2,:)    = cos(i*phi1-j*phi2);                
                br=br+2; 
        end
    end
%%

   
%evaluates the derivatives of base functions
function [v]=calculateV(phi1,phi2,K,bn,mr)
   v(K,length(phi1))=0; 
    
  br=1;  
  
  if mr==1 %partial derivatives is respect of the first variable
       for i=1:bn
                    v(br+1,:)    =  i*cos(i*phi1);
                    v(br+2,:)    = -i*sin(i*phi1);
                    br=br+2;
        end

        for i=1:bn
                    v(br+1,:)    = 0;
                    v(br+2,:)    = 0;
                    br=br+2;                
        end

        for i=1:bn
            for j=1:bn       
                    v(br+1,:)    = i*cos(i*phi1+j*phi2);               
                    v(br+2,:)    =-i*sin(i*phi1+j*phi2);
                    br=br+2; 
                    v(br+1,:)    = i*cos(i*phi1-j*phi2);                
                    v(br+2,:)    =-i*sin(i*phi1-j*phi2);                
                    br=br+2; 
            end
        end
    else %mr==2 and partial derivatives is respect of the second variable
        for i=1:bn
                    v(br+1,:)    = 0;
                    v(br+2,:)    = 0;
                    br=br+2;
        end

        for i=1:bn
                    v(br+1,:)    = i*cos(i*phi2);
                    v(br+2,:)    =-i*sin(i*phi2);
                    br=br+2;                
        end

        for i=1:bn
            for j=1:bn       
                    v(br+1,:)    = j*cos(i*phi1+j*phi2);               
                    v(br+2,:)    =-j*sin(i*phi1+j*phi2);
                    br=br+2; 
                    v(br+1,:)    =-j*cos(i*phi1-j*phi2);                
                    v(br+2,:)    = j*sin(i*phi1-j*phi2);                
                    br=br+2; 
            end
        end
    end
  
  
   
%%              

%evaluates the noise matrix
function [E]=calculateE(c,phiT,L,h,p)
    E=zeros(L,L);
    
    E=E+(phiT-c*p)*(phiT-c*p)';  
    E=(h/length(phiT(1,:)))*E;
%%   
  

%final result - evaluates the parameters
function [Cpt,XIpt]=calculateC(E,p,v1,v2,Cpr,XIpr,M,L,phiT,h)
    K=M/L;
    invr=inv(E);
    
    wr=lastwarn;
    if (~isempty(wr))&&((strcmp(wr(1:18),'Matrix is singular'))||(strcmp(wr(1:27),'Matrix is close to singular')))
        display('Singular matrix can lead to wrong and imprecise results. Please check the input parameters, signals or base functions. '); 
        error('Singular matrix can lead to wrong and imprecise results. Please check the input parameters, signals or base functions.');
    end
    
    
    %evaluates the concentration matrix ------
    XIpt=zeros(M,M);
        
    XIpt(1:K,1:K)        =XIpr(1:K,1:K)+h*invr(1,1)*p*(p');
    XIpt(1:K,K+1:2*K)    =XIpr(1:K,K+1:2*K)+h*invr(1,2)*p*(p');
    XIpt(K+1:2*K,1:K)    =XIpr(K+1:2*K,1:K)+h*invr(2,1)*p*(p');
    XIpt(K+1:2*K,K+1:2*K)=XIpr(K+1:2*K,K+1:2*K)+h*invr(2,2)*p*(p');
    
    
    %evaluates temp r ------
    r=zeros(K,L);
    ED=(E\phiT); %same as: ED=(inv(E)*phiT); (but faster)

    r(:,1)=XIpr(1:K,1:K)*Cpr(:,1)+XIpr(1:K,K+1:2*K)*Cpr(:,2)+h*(p*(ED(1,:)')-(1/2)*sum(v1,2));
    r(:,2)=XIpr(K+1:2*K,1:K)*Cpr(:,1)+XIpr(K+1:2*K,K+1:2*K)*Cpr(:,2)+h*(p*(ED(2,:)')-(1/2)*sum(v2,2));


%final evaluation of parameters 
 C=(XIpt\vertcat(r(:,1),r(:,2)))'; 
Cpt(:,1)=C(1:K); Cpt(:,2)=C(K+1:2*K);
%%

