alpha=0.6;
flag=0;
p=ones(n_a,n_g)/(15*pi);
p1=ones(n_a,n_g)/(15*pi);
pp=ones(n_a,n_g)/(15*pi);
pp1=ones(n_a,n_g)/(15*pi);
ppp=ones(n_a,n_g)/(15*pi);
ppp1=ones(n_a,n_g)/(15*pi);
% p(1:2,:)=0.016423994352368;

p(n_a-1:n_a,:)=0;  
p(:,1:2)=0;
p(:,n_g-1:n_g)=0;
p1=p;
pp=p;
pp1=p;
ppp=p;
ppp1=p;
for k=1:180000;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   p1(1,:)=0.436100618772761/6/pi;
   p1(2,:)=0.787519714482362/6/pi;
    for i=3:n_a-2;

        
        p1(i,1)=p(i,1); %p(i,n_g-3);% 
        p1(i,2)=p(i,2); % %p(i,n_g-2);
        for j=3:n_g-2;
            

% 
            c1=-s_a1(i-2,j)/12.0/delta_a^2-m_a1(i-2,j)/12.0/delta_a;
            c2=+16.0*s_a1(i-1,j)/12.0/delta_a^2+8.0*m_a1(i-1,j)/12.0/delta_a;
            c3=-s_g1(i,j-2)/12.0/delta_g^2-m_g1(i,j-2)/12.0/delta_g;
            c4=+16.0*s_g1(i,j-1)/12.0/delta_g^2+8.0*m_g1(i,j-1)/12.0/delta_g;
            c5=-30.0*s_a1(i,j)/12.0/delta_a^2-30.0*s_g1(i,j)/12.0/delta_g^2-1.4;%+beta/2.0+d/2.0/w^2/a(i)^2+e/2.0/w/a(i)*sin(gama(j));
            c6=16.0*s_g1(i,j+1)/12.0/delta_g^2-8.0*m_g1(i,j+1)/12.0/delta_g;
            c7=-s_g1(i,j+2)/12.0/delta_g^2+m_g1(i,j+2)/12.0/delta_g;
            c8=16.0*s_a1(i+1,j)/12.0/delta_a^2-8.0*m_a1(i+1,j)/12.0/delta_a;
            c9=-s_a1(i+2,j)/12.0/delta_a^2+m_a1(i+2,j)/12.0/delta_a;
            
            p1(i,j)=p(i,j)-alpha*(c1*p1(i-2,j)+c2*p1(i-1,j)+c3*p1(i,j-2)+c4*p1(i,j-1)+c5*p(i,j)+...  %
                c6*p(i,j+1)+c7*p(i,j+2)+c8*p(i+1,j)+c9*p(i+2,j)+1.50*pp(i,j)+1.50*ppp(i,j))/c5; %%%



            if(p1(i,j)<0.0)
                p1(i,j)=0.0;
            end
        end
     
      
        p1(i,n_g)=p1(i,4); 
        p1(i,n_g-1)=p1(i,3);
    end 
    p1(n_a-1:n_a,:)=0;
    
 %**********************************************
 
     
    pp1(1,:)=0.436100618772761/6/pi;;%0.0112974225120015/2/pi;
    pp1(2,:)=0.787519714482362/6/pi;%0.0112974225120015/2/pi;
    for i=3:n_a-2;

        
        pp1(i,1)=pp(i,1); %   %pp(i,n_g-3);
        pp1(i,2)=pp(i,2);% %pp(i,n_g-2);
        for j=3:n_g-2;
            

% 
            cc1=-s_a2(i-2,j)/12.0/delta_a^2-m_a2(i-2,j)/12.0/delta_a;
            cc2=+16.0*s_a2(i-1,j)/12.0/delta_a^2+8.0*m_a2(i-1,j)/12.0/delta_a;
            cc3=-s_g2(i,j-2)/12.0/delta_g^2-m_g2(i,j-2)/12.0/delta_g;
            cc4=+16.0*s_g2(i,j-1)/12.0/delta_g^2+8.0*m_g2(i,j-1)/12.0/delta_g;
            cc5=-30.0*s_a2(i,j)/12.0/delta_a^2-30.0*s_g2(i,j)/12.0/delta_g^2-3.0;%+beta/2.0+d/2.0/w^2/a(i)^2+e/2.0/w/a(i)*sin(gama(j));
            cc6=16.0*s_g2(i,j+1)/12.0/delta_g^2-8.0*m_g2(i,j+1)/12.0/delta_g;
            cc7=-s_g2(i,j+2)/12.0/delta_g^2+m_g2(i,j+2)/12.0/delta_g;
            cc8=16.0*s_a2(i+1,j)/12.0/delta_a^2-8.0*m_a2(i+1,j)/12.0/delta_a;
            cc9=-s_a2(i+2,j)/12.0/delta_a^2+m_a2(i+2,j)/12.0/delta_a;
            
            pp1(i,j)=pp(i,j)-alpha*(cc1*pp1(i-2,j)+cc2*pp1(i-1,j)+cc3*pp1(i,j-2)+cc4*pp1(i,j-1)+cc5*pp(i,j)+...  %
                cc6*pp(i,j+1)+cc7*pp(i,j+2)+cc8*pp(i+1,j)+cc9*pp(i+2,j)+0.7*p1(i,j)+1.5*ppp(i,j))/cc5; %


            if(pp1(i,j)<0.0)
                pp1(i,j)=0.0;
            end
        end
     
      
        pp1(i,n_g)=pp1(i,4);
        pp1(i,n_g-1)=pp1(i,3);
    end 
    pp1(n_a-1:n_a,:)=0;
    
    ppp1(1,:)=0.436100618772761/6/pi; %0.0112974225120015/2/pi;
    ppp1(2,:)=0.787519714482362/6/pi;%0.0112974225120015/2/pi;
    for i=3:n_a-2;

        
        ppp1(i,1)=ppp(i,1); %ppp(i,n_g-3); %
        ppp1(i,2)=ppp(i,2); %ppp(i,n_g-2); %
        for j=3:n_g-2;
            

% 
            cc1=-s_a3(i-2,j)/12.0/delta_a^2-m_a3(i-2,j)/12.0/delta_a;
            cc2=+16.0*s_a3(i-1,j)/12.0/delta_a^2+8.0*m_a3(i-1,j)/12.0/delta_a;
            cc3=-s_g3(i,j-2)/12.0/delta_g^2-m_g3(i,j-2)/12.0/delta_g;
            cc4=+16.0*s_g3(i,j-1)/12.0/delta_g^2+8.0*m_g3(i,j-1)/12.0/delta_g;
            cc5=-30.0*s_a3(i,j)/12.0/delta_a^2-30.0*s_g3(i,j)/12.0/delta_g^2-3.0;%+beta/2.0+d/2.0/w^2/a(i)^2+e/2.0/w/a(i)*sin(gama(j));
            cc6=16.0*s_g3(i,j+1)/12.0/delta_g^2-8.0*m_g3(i,j+1)/12.0/delta_g;
            cc7=-s_g3(i,j+2)/12.0/delta_g^2+m_g3(i,j+2)/12.0/delta_g;
            cc8=16.0*s_a3(i+1,j)/12.0/delta_a^2-8.0*m_a3(i+1,j)/12.0/delta_a;
            cc9=-s_a3(i+2,j)/12.0/delta_a^2+m_a3(i+2,j)/12.0/delta_a;
            
            ppp1(i,j)=ppp(i,j)-alpha*(cc1*ppp1(i-2,j)+cc2*ppp1(i-1,j)+cc3*ppp1(i,j-2)+cc4*ppp1(i,j-1)+cc5*ppp(i,j)+...  %
                cc6*ppp(i,j+1)+cc7*ppp(i,j+2)+cc8*ppp(i+1,j)+cc9*ppp(i+2,j)+0.7*p1(i,j)+1.5*pp1(i,j))/cc5;%

            if(ppp1(i,j)<0.0)
                ppp1(i,j)=0.0;
            end
        end
     
      
        ppp1(i,n_g)=ppp1(i,4);
        ppp1(i,n_g-1)=ppp1(i,3);
    end 
%%%%%%%归一化%%%%%%%%%%%%%%

if(mod(k,10)==0)
    s=0.0;
    for i=1:n_a
        for j=1:n_g
            s=s+(p1(i,j)+pp1(i,j)+ppp1(i,j))*delta_a*delta_g;
        end
    end
p1=p1./s;
pp1=pp1./s;
ppp1=ppp1./s;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=p1;
pp=pp1;
ppp=ppp1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_a=zeros(1,n_a);
for i=1:n_a
    for j=1:n_g
    p_a(1,i)=p_a(1,i)+(p(i,j)+pp(i,j)+ppp(i,j))*delta_g;
    end
end
figure
a=a;
plot(a,p_a)
hold on 
pa_moni=load('D:\c++\用户目录\我的文档\Visual Studio 2013\Projects\Pa_moni_paper1\Pa_moni_paper1\puth1ha1.2c0.008.txt');
plot(pa_moni(1:10:end,1),pa_moni(1:10:end,2),'o')
xlim([0 2.5])
hold off
%%%%%%%%%%%%%%%%%%%%%%
% n_x1=1000;
% n_x2=1000;
% x1=linspace(-2.0,2.0,n_x1);
% x2=linspace(-2.0,2.0,n_x2);
% delta_x1=x1(2)-x1(1);
% delta_x2=x2(2)-x2(1);
% p_xy=zeros(n_x1,n_x2);
% s=0.0;
% nn=0;
% for i=1:n_x1
%     for j=1:n_x2
%         s=x1(i).^2.0/2.0+0.3*x1(i).^4/4.0+x2(j).^2.0/2.0;
%         nn=floor(s/delta_a+1.0);
%         if(nn>n_a)
%          p_xy(i,j)=0.0;
%         else
%         p_xy(i,j)=p_a(nn)/t_h(nn,1);
%         end
%     end
% end
% s=0.0;
% for i=1:n_x1
%     for j=1:n_x2
%         s=s+p_xy(i,j)*delta_x1*delta_x2;
%     end
% end
% p_xy=p_xy./s;
% %%%%%%%%%%%%%%%%%%
% p_x=zeros(1,n_x1);
% for i=1:n_x1
%     for j=1:n_x2
%         p_x(1,i)=p_x(1,i)+p_xy(i,j)*delta_x2;
%     end
% end
% figure
% plot(x1,p_x);
% hold on
% px_moni=load('D:\c++\用户目录\我的文档\Visual Studio 2013\Projects\Pa_moni_paper1\Pa_moni_paper1\Px.txt');
% plot(px_moni(1:10:end,1),px_moni(1:10:end,2),'o')
% xlim([-2.5 2.5])
% hold off