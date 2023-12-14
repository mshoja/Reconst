function ModelPDF=FP(L,R,Mesh,D1,D2,DD2,PDF0,tau)
                
m = 0;
t = linspace(0,tau,3);
x=linspace(L,R,length(Mesh));
options = odeset('Stats','on','RelTol',1e-3,'AbsTol',1e-4);
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t,options);
% sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
u = sol(:,:,1);
ModelPDF=u(end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = pdex1pde(x,t,u,DuDx)

c = 1;

ndx=find(Mesh>=x,1);
d1=D1(ndx-1)+(D1(ndx)-D1(ndx-1))/(Mesh(ndx)-Mesh(ndx-1))*(x-Mesh(ndx-1));
d2=D2(ndx-1)+(D2(ndx)-D2(ndx-1))/(Mesh(ndx)-Mesh(ndx-1))*(x-Mesh(ndx-1));
dd2=DD2(ndx-1)+(DD2(ndx)-DD2(ndx-1))/(Mesh(ndx)-Mesh(ndx-1))*(x-Mesh(ndx-1));
 
f=-d1*u+dd2*u+d2*DuDx;
s = 0;
end
% ------------------------------------------------------------
function u0 = pdex1ic(x)
u0=PDF0(x);
end
% ------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
end
end
