clear;
clc;
format long;
n_a=250;n_g=300;
a=linspace(0.0000001,2.5,n_a);
gama=linspace(0,2*pi,n_g-3);
delta_a=a(2)-a(1);
delta_g=gama(2)-gama(1);
gama=[gama gama(n_g-3)+delta_g gama(n_g-3)+2*delta_g gama(n_g-3)+3*delta_g];
t_h=zeros(n_a,n_g);
g_h=zeros(n_a,n_g);
m_a1=zeros(n_a,n_g);
m_g1=zeros(n_a,n_g);
s_a1=zeros(n_a,n_g);
s_g1=zeros(n_a,n_g);

m_a2=zeros(n_a,n_g);
m_g2=zeros(n_a,n_g);
s_a2=zeros(n_a,n_g);
s_g2=zeros(n_a,n_g);

m_a3=zeros(n_a,n_g);
m_g3=zeros(n_a,n_g);
s_a3=zeros(n_a,n_g);
s_g3=zeros(n_a,n_g);

w=1.0;
omiga=1.2;

h1=1.0;
beta1=0.1;
e1=0.2;


% h2=0.95;
% beta2=0.102;
% e2=0.195;

h2=0.8;
beta2=0.11;
e2=0.18;

h3=0.8;
beta3=0.11;
e3=0.18;


d=0.004;
alfa=0.3;

%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n_a;
    for j=1:n_g;
        
% *EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE        
       a_h=sqrt((sqrt(w^4+4*alfa*a(i))-w^2)/alfa);
       f1=@(x)1./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       t_h(i,j)=4.0*quad(f1,0.0000000000001,a_h-0.000000000001);
       
       f2=@(x)-beta1.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))+e1.*sin(gama(j)).*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))./sqrt(2.0.*a(i))+d.*h1.^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       m_a1(i,j)=4.0*quad(f2,0.0000000000001,a_h-0.000000000001)/t_h(i,j);
       f3=@(x)((w.^2.*x+alfa.*x.^3)./sqrt(2*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))-e1.*cos(gama(j)).*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./sqrt(2)./a(i).^(1.5))./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
        m_g1(i,j)=omiga-4.0*quad(f3,0.0000000000001,a_h-0.000000000001)/t_h(i,j);
       f4=@(x)2.*d.*h1.^2.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
        s_a1(i,j)=2.0*quad(f4,0.000000000001,a_h-0.00000000001)/t_h(i,j);
        f5=@(x)d.*h1.^2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./a(i).^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
        s_g1(i,j)=2.0*quad(f5,0.000000000001,a_h-0.000000000001)/t_h(i,j);
% % *EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE               
       f2=@(x)-beta2.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))+e2.*sin(gama(j)).*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))./sqrt(2.0.*a(i))+d.*h2.^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       m_a2(i,j)=4.0*quad(f2,0.000000000001,a_h-0.000000000001)/t_h(i,j);
       f3=@(x)((w.^2.*x+alfa.*x.^3)./sqrt(2*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))-e2.*cos(gama(j)).*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./sqrt(2)./a(i).^(1.5))./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       m_g2(i,j)=omiga-4.0*quad(f3,0.000000000001,a_h-0.000000000001)/t_h(i,j);
       f4=@(x)2.*d.*h2.^2.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       s_a2(i,j)=2.0*quad(f4,0.000000000001,a_h-0.000000000001)/t_h(i,j);
       f5=@(x)d.*h2.^2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./a(i).^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       s_g2(i,j)=2.0*quad(f5,0.000000000001,a_h-0.000000000001)/t_h(i,j);
 %%8888888888888888888888888888888888888888888888888888888888
       f2=@(x)-beta3.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))+e3.*sin(gama(j)).*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))./sqrt(2.0.*a(i))+d.*h3.^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       m_a3(i,j)=4.0*quad(f2,0.000000000001,a_h-0.000000000001)/t_h(i,j);
       f3=@(x)((w.^2.*x+alfa.*x.^3)./sqrt(2*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))-e3.*cos(gama(j)).*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./sqrt(2)./a(i).^(1.5))./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       m_g3(i,j)=omiga-4.0*quad(f3,0.000000000001,a_h-0.000000000001)/t_h(i,j);
       f4=@(x)2.*d.*h3.^2.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       s_a3(i,j)=2.0*quad(f4,0.000000000001,a_h-0.000000000001)/t_h(i,j);
       f5=@(x)d.*h3.^2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./a(i).^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
       s_g3(i,j)=2.0*quad(f5,0.000000000001,a_h-0.000000000001)/t_h(i,j);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% alpha=0.4;
% flag=0;
% p=ones(n_a,n_g)/(5*pi);
% p1=ones(n_a,n_g)/(5*pi);
% pp=ones(n_a,n_g)/(5*pi);
% pp1=ones(n_a,n_g)/(5*pi);
% ppp=ones(n_a,n_g)/(5*pi);
% ppp1=ones(n_a,n_g)/(5*pi);
% % p(1:2,:)=0.016423994352368;
% 
% p(n_a-1:n_a,:)=0;  
% p(:,1:2)=0;
% p(:,n_g-1:n_g)=0;
% p1=p;
% pp=p;
% pp1=p;
% ppp=p;
% ppp1=p;
% for k=1:4000;
%     
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%    p1(1,:)=0.464765972222222/6/pi;
%    p1(2,:)=0.840345138888889/6/pi;
%     for i=3:n_a-2;
% 
%         
%         p1(i,1)=p(i,n_g-3);
%         p1(i,2)=p(i,n_g-2);
%         for j=3:n_g-2;
%             
% 
% % 
%             c1=-s_a1(i-2,j)/12.0/delta_a^2-m_a1(i-2,j)/12.0/delta_a;
%             c2=+16.0*s_a1(i-1,j)/12.0/delta_a^2+8.0*m_a1(i-1,j)/12.0/delta_a;
%             c3=-s_g1(i,j-2)/12.0/delta_g^2-m_g1(i,j-2)/12.0/delta_g;
%             c4=+16.0*s_g1(i,j-1)/12.0/delta_g^2+8.0*m_g1(i,j-1)/12.0/delta_g;
%             c5=-30.0*s_a1(i,j)/12.0/delta_a^2-30.0*s_g1(i,j)/12.0/delta_g^2-1.4;%+beta/2.0+d/2.0/w^2/a(i)^2+e/2.0/w/a(i)*sin(gama(j));
%             c6=16.0*s_g1(i,j+1)/12.0/delta_g^2-8.0*m_g1(i,j+1)/12.0/delta_g;
%             c7=-s_g1(i,j+2)/12.0/delta_g^2+m_g1(i,j+2)/12.0/delta_g;
%             c8=16.0*s_a1(i+1,j)/12.0/delta_a^2-8.0*m_a1(i+1,j)/12.0/delta_a;
%             c9=-s_a1(i+2,j)/12.0/delta_a^2+m_a1(i+2,j)/12.0/delta_a;
%             
%             p1(i,j)=p(i,j)-alpha*(c1*p1(i-2,j)+c2*p1(i-1,j)+c3*p1(i,j-2)+c4*p1(i,j-1)+c5*p(i,j)+...  %
%                 c6*p(i,j+1)+c7*p(i,j+2)+c8*p(i+1,j)+c9*p(i+2,j))/c5; %%%+1.50*pp(i,j)+1.50*ppp(i,j)
% 
% 
% 
%             if(p1(i,j)<0.0)
%                 p1(i,j)=0.0;
%             end
%         end
%      
%       
%         p1(i,n_g)=0; 
%         p1(i,n_g-1)=0;
%     end 
%     p1(n_a-1:n_a,:)=0;
%     
%  %**********************************************
%  
%      
%     pp1(1,:)=0.464765972222222/6/pi;
%     pp1(2,:)=0.840345138888889/6/pi;
%     for i=3:n_a-2;
% 
%         
%         pp1(i,1)=pp(i,n_g-3);
%         pp1(i,2)=pp(i,n_g-2);
%         for j=3:n_g-2;
%             
% 
% % 
%             cc1=-s_a2(i-2,j)/12.0/delta_a^2-m_a2(i-2,j)/12.0/delta_a;
%             cc2=+16.0*s_a2(i-1,j)/12.0/delta_a^2+8.0*m_a2(i-1,j)/12.0/delta_a;
%             cc3=-s_g2(i,j-2)/12.0/delta_g^2-m_g2(i,j-2)/12.0/delta_g;
%             cc4=+16.0*s_g2(i,j-1)/12.0/delta_g^2+8.0*m_g2(i,j-1)/12.0/delta_g;
%             cc5=-30.0*s_a2(i,j)/12.0/delta_a^2-30.0*s_g2(i,j)/12.0/delta_g^2-3.0;%+beta/2.0+d/2.0/w^2/a(i)^2+e/2.0/w/a(i)*sin(gama(j));
%             cc6=16.0*s_g2(i,j+1)/12.0/delta_g^2-8.0*m_g2(i,j+1)/12.0/delta_g;
%             cc7=-s_g2(i,j+2)/12.0/delta_g^2+m_g2(i,j+2)/12.0/delta_g;
%             cc8=16.0*s_a2(i+1,j)/12.0/delta_a^2-8.0*m_a2(i+1,j)/12.0/delta_a;
%             cc9=-s_a2(i+2,j)/12.0/delta_a^2+m_a2(i+2,j)/12.0/delta_a;
%             
%             pp1(i,j)=pp(i,j)-alpha*(cc1*pp1(i-2,j)+cc2*pp1(i-1,j)+cc3*pp1(i,j-2)+cc4*pp1(i,j-1)+cc5*pp(i,j)+...  %
%                 cc6*pp(i,j+1)+cc7*pp(i,j+2)+cc8*pp(i+1,j)+cc9*pp(i+2,j))/cc5; %+0.7*p1(i,j)+1.5*ppp(i,j)
% 
% 
%             if(pp1(i,j)<0.0)
%                 pp1(i,j)=0.0;
%             end
%         end
%      
%       
%         pp1(i,n_g)=0;
%         pp1(i,n_g-1)=0;
%     end 
%     pp1(n_a-1:n_a,:)=0;
%     
%     ppp1(1,:)=0.464765972222222/6/pi;
%     ppp1(2,:)=0.840345138888889/6/pi;
%     for i=3:n_a-2;
% 
%         
%         ppp1(i,1)=pp(i,n_g-3);
%         ppp1(i,2)=pp(i,n_g-2);
%         for j=3:n_g-2;
%             
% 
% % 
%             cc1=-s_a3(i-2,j)/12.0/delta_a^2-m_a3(i-2,j)/12.0/delta_a;
%             cc2=+16.0*s_a3(i-1,j)/12.0/delta_a^2+8.0*m_a3(i-1,j)/12.0/delta_a;
%             cc3=-s_g3(i,j-2)/12.0/delta_g^2-m_g3(i,j-2)/12.0/delta_g;
%             cc4=+16.0*s_g3(i,j-1)/12.0/delta_g^2+8.0*m_g3(i,j-1)/12.0/delta_g;
%             cc5=-30.0*s_a3(i,j)/12.0/delta_a^2-30.0*s_g3(i,j)/12.0/delta_g^2-3.0;%+beta/2.0+d/2.0/w^2/a(i)^2+e/2.0/w/a(i)*sin(gama(j));
%             cc6=16.0*s_g3(i,j+1)/12.0/delta_g^2-8.0*m_g3(i,j+1)/12.0/delta_g;
%             cc7=-s_g3(i,j+2)/12.0/delta_g^2+m_g3(i,j+2)/12.0/delta_g;
%             cc8=16.0*s_a3(i+1,j)/12.0/delta_a^2-8.0*m_a3(i+1,j)/12.0/delta_a;
%             cc9=-s_a3(i+2,j)/12.0/delta_a^2+m_a3(i+2,j)/12.0/delta_a;
%             
%             ppp1(i,j)=ppp(i,j)-alpha*(cc1*ppp1(i-2,j)+cc2*ppp1(i-1,j)+cc3*ppp1(i,j-2)+cc4*ppp1(i,j-1)+cc5*ppp(i,j)+...  %
%                 cc6*ppp(i,j+1)+cc7*ppp(i,j+2)+cc8*ppp(i+1,j)+cc9*ppp(i+2,j))/cc5;%+0.7*p1(i,j)+1.5*pp1(i,j)
% 
%             if(ppp1(i,j)<0.0)
%                 ppp1(i,j)=0.0;
%             end
%         end
%      
%       
%         ppp1(i,n_g)=0;
%         ppp1(i,n_g-1)=0;
%     end 
% %%%%%%%%¹éÒ»»¯%%%%%%%%%%%%%%
% 
% if(mod(k,10)==0)
%     s=0.0;
%     for i=1:n_a
%         for j=1:n_g
%             s=s+(p1(i,j)+pp1(i,j)+ppp1(i,j))*delta_a*delta_g;
%         end
%     end
% p1=p1./s;
% pp1=pp1./s;
% ppp1=ppp1./s;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p=p1;
% pp=pp1;
% ppp=ppp1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_a=zeros(1,n_a);
% for i=1:n_a
%     for j=1:n_g
%     p_a(1,i)=p_a(1,i)+(p(i,j)+pp(i,j)+ppp(i,j))*delta_g;
%     end
% end
% 
% %clc;
% % format long;
% % n_a=200;n_g=160;
% % a=linspace(0.000001,3,n_a);
% % gama=linspace(0,2*pi,n_g-3);
% % delta_a=a(2)-a(1);
% % delta_g=gama(2)-gama(1);
% % gama=[gama gama(n_g-3)+delta_g gama(n_g-3)+2*delta_g gama(n_g-3)+3*delta_g];
% % t_h=zeros(n_a,n_g);
% % g_h=zeros(n_a,n_g);
% % m_a1=zeros(n_a,n_g);
% % m_g1=zeros(n_a,n_g);
% % s_a1=zeros(n_a,n_g);
% % s_g1=zeros(n_a,n_g);
% % 
% % m_a2=zeros(n_a,n_g);
% % m_g2=zeros(n_a,n_g);
% % s_a2=zeros(n_a,n_g);
% % s_g2=zeros(n_a,n_g);
% % 
% % m_a3=zeros(n_a,n_g);
% % m_g3=zeros(n_a,n_g);
% % s_a3=zeros(n_a,n_g);
% % s_g3=zeros(n_a,n_g);
% % 
% % w=1.0;
% % omiga=1.2;
% % h1=1.0;
% % beta1=0.1;
% % e1=0.2;
% % h2=0.95;
% % beta2=0.102;
% % e2=0.195;
% % h3=0.8;
% % beta3=0.11;
% % e3=0.18;
% % 
% % d=0.004;
% % alfa=0.3;
% % 
% % %%%%%%%%%%%%%%
% % % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%
% % for i=1:n_a;
% %     for j=1:n_g;
% %         
% % % *EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE        
% %        a_h=sqrt((sqrt(w^4+4*alfa*a(i))-w^2)/alfa);
% %        f1=@(x)1./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %        t_h(i,j)=4.0*quad(f1,0.0000000000001,a_h-0.000000000001);
% %        
% %        f2=@(x)-beta1.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))+e1.*sin(gama(j)).*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))./sqrt(2.0.*a(i))+d.*h1.^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %        m_a1(i,j)=4.0*quad(f2,0.0000000000001,a_h-0.000000000001)/t_h(i,j);
% %        f3=@(x)((w.^2.*x+alfa.*x.^3)./sqrt(2*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))-e1.*cos(gama(j)).*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./sqrt(2)./a(i).^(1.5))./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %         m_g1(i,j)=omiga-4.0*quad(f3,0.0000000000001,a_h-0.000000000001)/t_h(i,j);
% %        f4=@(x)2.*d.*h1.^2.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %         s_a1(i,j)=2.0*quad(f4,0.000000000001,a_h-0.00000000001)/t_h(i,j);
% %         f5=@(x)d.*h1.^2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./a(i).^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %         s_g1(i,j)=2.0*quad(f5,0.000000000001,a_h-0.000000000001)/t_h(i,j);
% % % % *EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE               
% %        f2=@(x)-beta2.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))+e2.*sin(gama(j)).*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))./sqrt(2.0.*a(i))+d.*h2.^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %        m_a2(i,j)=4.0*quad(f2,0.000000000001,a_h-0.000000000001)/t_h(i,j);
% %        f3=@(x)((w.^2.*x+alfa.*x.^3)./sqrt(2*(w.^2./2.0*x.^2+alfa./4.0.*x.^4))-e2.*cos(gama(j)).*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./sqrt(2)./a(i).^(1.5))./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %        m_g2(i,j)=omiga-4.0*quad(f3,0.000000000001,a_h-0.000000000001)/t_h(i,j);
% %        f4=@(x)2.*d.*h2.^2.*sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %        s_a2(i,j)=2.0*quad(f4,0.000000000001,a_h-0.000000000001)/t_h(i,j);
% %        f5=@(x)d.*h2.^2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4)./a(i).^2./sqrt(2.*a(i)-2.*(w.^2./2.0*x.^2+alfa./4.0.*x.^4));
% %        s_g2(i,j)=2.0*quad(f5,0.000000000001,a_h-0.000000000001)/t_h(i,j);
% %  %%8888888888888888888888888888888888888888888888888888888888
% %        f