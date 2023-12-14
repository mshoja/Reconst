function x=SimMayLevy(x0,dt,T,r,K,AL,GA,alpha,betta,sigma)

x=zeros(1,T,'single');x(1)= x0;
pd=makedist('Stable','alpha',alpha,'beta',betta,'gam',1/sqrt(2),'delta',0);Noise=single(dt^(1/alpha)*(sqrt(2*sigma))*random(pd,[1 T]));
clear pd

for i=2:T
    % Typical integration
    x(i)=x(i-1)+(r*x(i-1)*(1-x(i-1)/K)-GA*x(i-1)^2/(AL^2+x(i-1)^2))*dt+Noise(i);

    % Heuns integration method (more accurate)
    % z1=r*x(i-1)*(1-x(i-1)/K)-GA*x(i-1)^2/(AL^2+x(i-1)^2);
    % z2=x(i-1)+z1*dt;
    % x(i)=x(i-1)+(z1+r*z2*(1-z2/K)-GA*z2^2/(AL^2+z2^2))*dt/2+Noise(i);

    % if x(i)<=0, break;end
end

clear Noise
end
