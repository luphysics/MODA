function [t1,t2,q1,q2]=CFprint(cc,bn)
%plots the coupling functions from the inferred parameters

%---inputs---
%cc - vector of inferred parameters
%bn - order of Fourier base function

%Note that the input is vector of parameters for one time window
%%
            %---evaluating the coupling functions -----
            t1=0:0.13:2*pi;t2=0:0.13:2*pi; 
            q1(1:length(t1),1:length(t1))=0;q2=q1;
            u=cc; K=length(u)/2;
            for i1=1:length(t1)                
                for j1=1:length(t2)
                    br=2;
                   
                    for ii=1:bn
                        q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t1(i1))+u(br+1)*cos(ii*t1(i1));
                        q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(ii*t2(j1))+u(K+br+1)*cos(ii*t2(j1));
                        br=br+2;  
                    end
                    for ii=1:bn
                        q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t2(j1))+u(br+1)*cos(ii*t2(j1));
                        q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(ii*t1(i1))+u(K+br+1)*cos(ii*t1(i1));
                        br=br+2;
                    end
                    
                    for ii=1:bn
                        for jj=1:bn                            
                                   q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t1(i1)+jj*t2(j1))+u(br+1)*cos(ii*t1(i1)+jj*t2(j1));                                                                
                                   q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(ii*t1(i1)+jj*t2(j1))+u(K+br+1)*cos(ii*t1(i1)+jj*t2(j1));                                   
                                   br=br+2;
                                   
                                   q1(i1,j1)=q1(i1,j1)+u(br)*sin(ii*t1(i1)-jj*t2(j1))+u(br+1)*cos(ii*t1(i1)-jj*t2(j1));                                                                
                                   q2(i1,j1)=q2(i1,j1)+u(K+br)*sin(ii*t1(i1)-jj*t2(j1))+u(K+br+1)*cos(ii*t1(i1)-jj*t2(j1));     
                                   br=br+2;
                        end
                    end                    
                                               
                end
            end

                        %---plotting -----
%                         f1=figure;
% 
%                         subplot(1,2,1);surf(t1,t2,q1','FaceColor','interp');                                              
%                         view([-40 50])
%                         set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')
%                         xlabel('\phi_1');ylabel('\phi_2');zlabel('q_1(\phi_1,\phi_2)');axis tight
% 
%                         subplot(1,2,2);surf(t1,t2,q2','FaceColor','interp');                                               
%                         view([-40 50])
%                         set(gca,'fontname','Helvetica','fontsize',12,'Xgrid','off','Ygrid','off')            
%                         xlabel('\phi_1');ylabel('\phi_2');zlabel('q_2(\phi_1,\phi_2)');axis tight
%                         
%                         colormap(hot)
%                         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

                         %uncomment this lines for saving the figure
                          % saveas(f1,'filename','jpg');
                          % saveas(f1,'filename','fig');

%% 
