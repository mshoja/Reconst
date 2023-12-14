function x=SimQuinticLevy(L,R,x0,dt,T,alpha,betta,sigma)

x=zeros(1,T,'single');x(1)= x0;
pd=makedist('Stable','alpha',alpha,'beta',betta,'gam',1/sqrt(2),'delta',0);Noise=single(random(pd,[1 T]));
clear pd

for i=2:T
    if x(i-1)>L && x(i-1)<R
%         x(i)=x(i-1)+(-x(i-1)^3)*dt+dt^(1/alpha)*(sqrt(2*sigma))*Noise(i);
%         % Heuns integration method 
%           z=x(i-1)+(-x(i-1)^3)*dt;
%           x(i)=x(i-1)+(-x(i-1)^3-z^3)*dt/2+dt^(1/alpha)*(sqrt(2*sigma))*Noise(i);
          %Heuns integration with Horner representation (more accurate)
          x(i)=-x(i-1)*(x(i-1)^2*(dt + x(i-1)^2*(x(i-1)^2*(- (x(i-1)^2*dt^4)/2 + (3*dt^3)/2) - (3*dt^2)/2)) - 1)+dt^(1/alpha)*(sqrt(2*sigma))*Noise(i); 
    elseif x(i-1)<L
        x(i)=L+dt;
    else
        x(i)=R-dt;
    end
end
clear Noise

