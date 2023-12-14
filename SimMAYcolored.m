
function x=SimMAYcolored(L,R,r,K,AL,GA,x0,dt,T,D,tau)

x = zeros(1,T,'single');xi=x;x(1)= x0; 
switch tau
    case 0
        NoiseW=single(sqrt(2*D*dt)*randn(1,T));
        for i=2:T
            if (x(i-1)>L && x(i-1)<R)
                x(i)=x(i-1)+(r*x(i-1)*(1-x(i-1)/K)-GA*x(i-1)^2/(AL^2+x(i-1)^2))*dt+NoiseW(i);
            elseif x(i-1)<L
                   x(i)=L+dt;
            else
               x(i)=R-dt;
            end
        end
    otherwise
        NoiseC=single(sqrt(2*dt*D)/tau*randn(1,T));
        for i=2:T
            if (x(i-1)>L && x(i-1)<R)
                xi(i)=xi(i-1)+(-xi(i-1)/tau)*dt+NoiseC(i);
                x(i)=x(i-1)+(r*x(i-1)*(1-x(i-1)/K)-GA*x(i-1)^2/(AL^2+x(i-1)^2)+xi(i-1))*dt;
            elseif x(i-1)<L
                x(i)=L+dt;
            else
                x(i)=R-dt;
            end
        end
end
clear xi NoiseW NoiseC   
    





