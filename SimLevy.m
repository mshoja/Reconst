function x=SimLevy(L,R,T,g,h,M,dt,x0,alpha,betta)

x=zeros(1,T,'single');x(1)= x0;
pd=makedist('Stable','alpha',alpha,'beta',betta,'gam',1,'delta',0);Noise=single(dt^(1/alpha)*random(pd,[1 T]));
clear pd

for i=2:T
    %Euler–Maruyama
    m=ceil((x(i-1)-L)./(R-L)*(M-1)+1);
    a=x(i-1)+g(m)*dt+h(m)*Noise(i); 
    a(a<=L)=L+dt;a(a>=R)=R-dt;
    x(i)=a;
    %Heuns integration method (more accurate though a bit more time-consumming)
%     m1=ceil((x(i-1)-L)./(R-L)*(M-1)+1);
%     z=x(i-1)+g(m1)*dt;
%     m2=ceil((z-L)./(R-L)*(M-1)+1);m2(m2<1)=1;m2(m2>M)=M;
%     a=x(i-1)+(g(m1)+g(m2))*dt/2+h(m1)*Noise(i);
%     a(a<=L)=L+dt;a(a>=R)=R-dt;
%     x(i)=a;
end
clear Noise
end


