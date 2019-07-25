function [ s ] = sync_map (c,bn)
% Calculates synchronization through map representation given the inferred parameters. 
% It applies to 1:1 synchonization ratio, higher n:m ratio should be
% adjust by deviding the phases by n and m beforehand. 
%The algorithmn employs modified Newton-Raphson root method. 


%---inputs---
%c  - vector of inferred parameters
%bn - order of Fourier base function

%---outputs---
%s - syncronization index: s=1 sync, s=0 nosync

%%
	s=0; 	
	eps=1.e-8;	h_rk=0.05;	
	max_step=0.2;	k_delta=0.1;	min_step=0.0001;

h_backward = -min_step;   
k_backward = -2*min_step; 

xk=0; 
xh=min_step; 
x0=2*min_step;

        fxk=unwrapped_map(xk,c,h_rk,eps,bn);

        fxh=unwrapped_map(xh,c,h_rk,eps,bn);

	if ((fxk-xk)*(fxh-xh)<0)  s = 1; return; end;
        fx0=unwrapped_map(x0,c,h_rk,eps,bn);

	if ((fxh-xh)*(fx0-x0)<0)  s = 1; return; end;
 	i=3;		

while(x0<2*pi+4*min_step)

	dx0 = uneven_forward_wrapping_aware_f2_differences(fx0,fxh,fxk,h_backward,k_backward);

	incre =   min( max( min_step,0.5*abs( wrapToPi(fx0 - x0)/ dx0(1) ) ) , max_step );
	
	k_backward = h_backward - incre;	
	h_backward = -incre;
	xk=xh;	fxk=fxh;
	xh=x0;	fxh=fx0;

	x0=x0+incre; fx0=unwrapped_map(x0,c,h_rk,eps,bn); 
	if ((fxh-xh)*(fx0-x0)<0)  s = 1; return; end;

end
end

%% -------------------------- unwrapped_map--------------------------------
function [ Dth_next ] = unwrapped_map (Dth,c,h,eps,bn);
th=[-Dth/2,Dth/2];


br=0;
br1=0;

th_next = rk4(th,c,h,bn);
th_old=th_next;
mth=th_next(1)+th_next(2);

while(abs(mth)<4*pi)&&(br1<500)
	th_old=th_next;
	th_next = rk4(th_next,c,h,bn);
	mth=th_next(1)+th_next(2);
    br1=br1+1;
end

if (br1==500) Dth_next=pi; return; end 


%  Newton-Raphson method
if (mth>0)
	while ((abs(mth)-4*pi)> eps)&&(br<100)
		dot_th_next= calculate_dot_th1th2(th_next,c,bn);
		h=  - (mth-4*pi)/(dot_th_next(1)+dot_th_next(2)); 
		th_next = rk4(th_next,c,h,bn);
		mth=th_next(1)+th_next(2);
        
         br=br+1;
	end
else
	while ((abs(mth)+4*pi)> eps)&&(br<100)
		dot_th_next= calculate_dot_th1th2(th_next,c,bn);
		h=  - (mth+4*pi)/(dot_th_next(1)+dot_th_next(2)); 
		th_next = rk4(th_next,c,h,bn);
		mth=th_next(1)+th_next(2);
        
         br=br+1;
	end

end

Dth_next=   th_next(2)-th_next(1)  ;

end

%% ------------------------------rk4---------------------------------------
function [ th_next ] = rk4(th_0,c,h,bn)

	k_1= calculate_dot_th1th2(th_0,c,bn);
	k_2= calculate_dot_th1th2((th_0 + 0.5*h*k_1),c,bn);
	k_3= calculate_dot_th1th2((th_0 + 0.5*h*k_2),c,bn);
	k_4= calculate_dot_th1th2((th_0 + h*k_3),c,bn);
	th_next=th_0 + h/6. * (k_1 + k_2 + k_3 + k_4);
end


%% -----------------------calculate_dot_th1th2-----------------------------
function [ dot_th1th2 ] = calculate_dot_th1th2(th,c,bn)
                   
                    K=length(c)/2;
                    dot_th1th2(1) = c(1);
                    dot_th1th2(2) = c(K+1);
                    br=2;
                   
                    for ii=1:bn
                        dot_th1th2(1)=dot_th1th2(1)+c(br)*sin(ii*th(1))+c(br+1)*cos(ii*th(1));
                        dot_th1th2(2)=dot_th1th2(2)+c(K+br)*sin(ii*th(2))+c(K+br+1)*cos(ii*th(2));
                        br=br+2;
                    end
                    
                    for ii=1:bn
                        dot_th1th2(1)=dot_th1th2(1)+c(br)*sin(ii*th(2))+c(br+1)*cos(ii*th(2));
                        dot_th1th2(2)=dot_th1th2(2)+c(K+br)*sin(ii*th(1))+c(K+br+1)*cos(ii*th(1));
                        br=br+2;
                    end
                    
                    for ii=1:bn
                        for jj=1:bn                            
                                   dot_th1th2(1)=dot_th1th2(1)+c(br)*sin(ii*th(1)+jj*th(2))+c(br+1)*cos(ii*th(1)+jj*th(2));                                                                
                                   dot_th1th2(2)=dot_th1th2(2)+c(br+1)*sin(ii*th(1)+jj*th(2))+c(br+3)*cos(ii*th(1)+jj*th(2));                                   
                                   br=br+2;
                                   
                                   dot_th1th2(1)=dot_th1th2(1)+c(br)*sin(ii*th(1)-jj*th(2))+c(br+1)*cos(ii*th(1)-jj*th(2));                                                                
                                   dot_th1th2(2)=dot_th1th2(2)+c(K+br)*sin(ii*th(1)-jj*th(2))+c(K+br+1)*cos(ii*th(1)-jj*th(2));    
                                   br=br+2;
                        end
                    end        

end

%% --------------uneven_forward_wrapping_aware_f2_differences--------------

function [ f ] = uneven_forward_wrapping_aware_f2_differences(fx0,fxh,fxk,h,k)

	fx0 = wrapToPi(fx0);
	fxh = wrapToPi(fxh);
	fxk = wrapToPi(fxk);
 
	temp= ( h*fxk  -k*fxh + (k-h)*fx0 ) / ( k*(k - h) );
	f(2)= 2*temp/h; 
	f(1)= (fxh - fx0)/h  - h*temp; 
end

%% -------------------------wrapToPi_--------------------------------------
function fx0 = wrapToPi(fx1)
fx0=mod(fx1,2*pi);
end
%% -------------------------END--------------------------------------------