
clc;

%NOTE: You need to download CHEBFUN from the following link (and then add it to your MATLAB working path)
%  https://nl.mathworks.com/matlabcentral/fileexchange/47023-chebfun-current-version 

%%
%*******************
%**** FIGURE 2A ****
%*******************

al=0.17;g=0.28;
% model parameters (for the model of May) 
r=0.1;K=20;h=al*K;ga=r*g*K;   % r is growth rate, K is carying capacity, h is half saturation constant, ga is consumption rate 

L=0;R=K;  % Lower and upper bounds considered for the state space
chebfunpref.setDefaults('domain',[L R]);
f = @(x) r*x*(1-x/K)-ga*x^2/(h^2+x^2);
f=chebfun(f,[L R],'vectorize');
E=roots(f);
u=chebfun(-cumsum(f),'vectorize');u=u-min(u);

D=[0.03 0.04 0.06 0.102]; % These are different noise intensities in Figure 2A. Note that noise-induced bifurcation occurs at D~0.102 
mesh=linspace(L,R,2000);   % A fine mesh over the state space
fig1=figure(1);
for d=D
    D1=r.*mesh.*(1-mesh./K)-ga.*mesh.^2./(h.^2+mesh.^2);
    D2=0.*D1+1./2.*(d.*mesh.*(1-mesh./K)).^2;DD2=0.*D1+d.^2.*mesh.*(1-2.*mesh./K).*(1-mesh./K);
    tf=50000;   % tf is the final time we use to solve Fokker-Plank equation to get stationar density (in theory it should be infinite) 
   
    si=0.1;
    PDF0=@(x)1/sqrt(2*pi*si^2)*exp(-0.5*((x-mean([L R]))/si)^2); % An initial Gaussian distribution placed at the center of state space and with std of si=0.1
    MayPDF=FP(L,R,mesh,D1,D2,DD2,PDF0,tf);  % FP numerically solves the Fokker-Planck equation (only for Gaussian white noise)

    switch d
        case D(1)
            plot(mesh,MayPDF,'-k');hold on;f1=chebfun(MayPDF',[L R],'equi',2000);
        case D(2)
            my_colour = [141 52 16] ./ 255;plot(mesh,MayPDF,'-','color',my_colour);hold on;
        case D(3)
            my_colour = [169 252 26] ./ 255;hold on;plot(mesh,MayPDF,'-','color','b');hold on;
        case D(4)
            my_colour = [35 148 63] ./ 255;hold on;plot(mesh,MayPDF,':','color',my_colour,'Linewidth',1.5);hold on;
    end
end
ylim([0 0.65]);xlim([0 16]);
xlabel('Ecosystem state (vegetation biomass)');ylabel('Stationary distribution');
title('Grazing model of May (bistable)');
hl1=line([E(2) E(2)],[0 f1(E(2))],'color','k','LineStyle','-.');hold on;hl2=line([E(4) E(4)],[0 0.063],'color','k','LineStyle','-.');hold on;
hasbehavior(hl1,'legend',false);hasbehavior(hl2,'legend',false);  

fig2=figure(2);
x2=linspace(0,K,1000);
y2=1./2.*(x2.*(1-x2./K)).^2;f2=chebfun(@(x)1/2*(x*(1-x/K))^2,[0 K],'vectorize');
hold on;

for d=D
    switch d
        case D(1)
            plot(d*f2,'-k');hold on;
        case D(2)
            my_colour = [141 52 16] ./ 255;plot(d*f2,'-','color',my_colour);hold on;
        case D(3)
            my_colour = [169 252 26] ./ 255;hold on;plot(d*f2,'-','color','b');hold on;
        case D(4)
            my_colour = [35 148 63] ./ 255;hold on;plot(d*f2,'-','color',my_colour,'LineStyle',':','Linewidth',1.5);hold on;
    end
end
f=D(4)*f2;
line([E(2) E(2)],[0 f(E(2))],'color','k','LineStyle','-.');hold on;hl2=line([E(4) E(4)],[0 f(E(4))],'color','k','LineStyle','-.');hold on;

box on;
set(gca,'fontsize',8);
title('Noise intensity');
[h_m h_i]=inset(fig1,fig2,0.35);
close(fig1,fig2);

%%
%*******************
%**** FIGURE 2B **** 
%*******************

% Note: This section is going to take around 7 min

al=0.17;g=0.28;
% model parameters (for the model of May) 
r=0.1;K=20;h=al*K;ga=r*g*K;   % r is growth rate, K is carying capacity, h is half saturation constant, ga is consumption rate 

L=0;R=14;  % Lower and upper bounds considered for the state space
chebfunpref.setDefaults('domain',[L R]);
f = @(x) r*x*(1-x/K)-ga*x^2/(h^2+x^2);
f=chebfun(f,[L R],'vectorize');
E=roots(f);
u=chebfun(-cumsum(f),'vectorize');u=u-min(u);

dt=10^(-2);T=10^7;x0=E(3);D=0.02;  % D is noise intensity  

Tau=[0 10 50 100 150 200 300 500];  % Tau is a vector of noise colours we use
DATA=zeros(2^12,8);i=1;
for tau=Tau
    y=single([]);
    parfor n=1:1500   % Note: if you choose a bigger n then the results are more accurate (at the expense of higher computational time). For the paper we used n=1500
        x=SimMAYcolored(L,R,r,K,h,ga,x0,dt,T,D,tau);y=[y x(round(0.05*T):100:end)];x=[];
    end
    [~,f,yi,~]=kde(y',2^12,L,R);
    DATA(:,i)=f;
    clear y
    i=i+1;
end

LW='LineWidth';
h=zeros(1,8);
for tau=Tau
    switch tau
        case Tau(1)
            fx=chebfun(DATA(:,1),[L R],'equi',2000);
            h(1)=plot(fx,'-k',LW,0.1);hold on;xlim([L R]);hold on;
        case Tau(2)
            my_colour = [141 52 16] ./ 255;
            fx=chebfun(DATA(:,2),[L R],'equi',2000);
            h(2)=plot(fx,LW,0.1,'color',my_colour);hold on;xlim([L R]);hold on;
        case Tau(3)
            my_colour = [35 148 63] ./ 255;
            fx=chebfun(DATA(:,3),[L R],'equi',2000);
            h(3)=plot(fx,LW,0.1,'color',my_colour);hold on;xlim([L R]);hold on;
        case Tau(4)
            fx=chebfun(DATA(:,4),[L R],'equi',2000);
            h(4)=plot(fx,LW,0.1,'color','b');hold on;xlim([L R]);hold on;
        case Tau(5)
            my_colour = [169 252 26] ./ 255;
            fx=chebfun(DATA(:,5),[L R],'equi',2000);
            h(5)=plot(fx,LW,0.1,'color',my_colour);hold on;xlim([L R]);hold on;
        case Tau(6)
            fx=chebfun(DATA(:,6),[L R],'equi',2000);
            h(6)=plot(fx,LW,0.1,'color','green');hold on;xlim([L R]);hold on;
        case Tau(7)
            fx=chebfun(DATA(:,7),[L R],'equi',2000);
            my_colour = [252 147 49] ./ 255;
            h(7)=plot(fx,LW,0.1,'color',my_colour);hold on;xlim([L R]);hold on;
        case Tau(8)
            fx=chebfun(DATA(:,8),[L R],'equi',2000);
            h(8)=plot(fx,'-r',LW,0.1);hold on;xlim([L R]);hold on;
    end
end
hold on;xlabel('Ecosystem state (vegetation biomass)');hold on;ylabel('Stationary density');hold on;ylim([0 0.7]);hold on;
title('Grazing model of May (bistable)');
line([E(2) E(2)],[0 max(fx)],'color','k','LineStyle','-.');hold on;line([E(4) E(4)],[0 max(fx)],'color','k','LineStyle','-.');
legend(h,['noise color =' num2str(0) ' (White noise)'],['noise color = ' num2str(Tau(2))],['noise color = ' num2str(Tau(3))],['noise color = ' num2str(Tau(4))],['noise color = ' num2str(Tau(5))],['noise color = ' num2str(Tau(6))],['noise color = ' num2str(Tau(7))],['noise color = ' num2str(Tau(8))]);

%%
%*******************
%**** FIGURE 2C **** 
%*******************
% Note: This section is going to take around 14 min

al=0.17;g=0.28;
% model parameters (for the model of May) 
r=0.1;K=20;h=al*K;ga=r*g*K;   % r is growth rate, K is carying capacity, h is half saturation constant, ga is consumption rate 

chebfunpref.setDefaults('factory');
L=0;R=100;
chebfunpref.setDefaults('domain',[L R]);
f = @(x) r*x*(1-x/K)-ga*x^2/(h^2+x^2);f=chebfun(f,'vectorize');E=roots(f);
u=chebfun(-cumsum(f),'vectorize');u=u-min(u);

dt=10^(-2);T=10^6;x0=E(3);
betta=0.03;sigma=0.01;  %( For other choice of parameters betta=0.02)
alpha=[1.95 1.85 1.75];


% Note: The following loop is memory conssuming. If your run out of memory consider reducing 'T'
% but to keep the same level of accuracy increase 'i' proportionally (For instance consider T=T/2, i=2*i) 

DATA=cell(1,3);
for n=1:3
    y=[];y1=SimMayLevy(x0,dt,T,r,K,h,ga,alpha(n),betta,sigma);
    for i=2:1500   % Note: to get a highly accurate result we need a big range for 'i'.   
        m=find(y1<=0,1);
        if isempty(m)
            y=[y y1(round(0.1*T):100:end)];   % round(0.1*T) is used to remove transient effects. WE take every 100 of data to reduce the computational burden of estimating the density
            y1=SimMayLevy(y(end),dt,T,r,K,h,ga,alpha(n),betta,sigma);
        else
            y=[y y1(round(0.1*T):100:m-1)];
            y1=SimMayLevy(x0,dt,T,r,K,h,ga,alpha(n),betta,sigma);
        end
    end
    clear y1
    DATA{n}=y;
end

a=[];
mesh=linspace(0,20,3000);
H=zeros(1,4);
for n=1:3
    y=DATA{n};
    y=y(y>0);y=y(y<20);
    [f,~] = ksdensity(y(1:end),mesh);
    a=max([a f]);
    switch n
        case 1
            my_colour = [141 52 16] ./ 255;H(2)=plot(mesh,f,'color',my_colour);hold on;
        case 2
            H(3)=plot(mesh,f,'-b');hold on;
        case 3
            my_colour = [35 148 63] ./ 255;H(4)=plot(mesh,f,'color',my_colour,'LineStyle',':','LineWidth',1.5);hold on;
    end
end

D1=r.*mesh.*(1-mesh./K)-ga.*mesh.^2./(h.^2+mesh.^2);
D2=0.*D1+sigma;DD2=0.*D1;
tf=50000;   % tf is the final time we use to solve Fokker-Plank equation to get stationar density (in theory it should be infinite) 
si=0.1;
PDF0=@(x)1./sqrt(2.*pi.*si.^2).*exp(-0.5.*((x-mean([0 20]))./si).^2); % An initial Gaussian distribution centered at the state space and with std of si=0.1
MayPDF=FP(0,20,mesh,D1,D2,DD2,PDF0,tf);  % FP numerically solves the Fokker-Planck equation (only for Gaussian white noise)
H(1)=plot(mesh,MayPDF,'-k');
xlabel('Ecosystem state (vegetation biomass)');ylabel('Stationary distribution');title('Grazing model of May (bistable)');
b=max(MayPDF);
line([E(2) E(2)],[0 b],'LineStyle','-.','color','k');line([E(4) E(4)],[0 a],'LineStyle','-.','color','k');hold on;
xlim([0 16]);ylim([0 0.4]);
legend(H,{'\alpha = 2 (Gaussian noise)','\alpha = 1.95','\alpha = 1.85','\alpha = 1.75'});
%%
%*******************
%**** FIGURE 2D **** 
%*******************

% Note: This section is going to take around 1 min

%I)Quartic potential 
clc;
chebfunpref.setDefaults('factory');
L=-40;R=40;chebfunpref.setDefaults('domain',[L R]);
alpha=1;betta=0;mu=0;sigma=0.5;x0=0;dt=10^(-3);T =2*10^(4); 

y=single([]);
parfor n=1:100000   % A bigger range for n leads to more accurate results at the expense of higher computational time (for the paper we used 200000)
    x=SimQuinticLevy(L,R,x0,dt,T,alpha,betta,sigma);
    y=[y x(round(0.1*T):100:end)];
    x=[];
end
a=200;
y=y(y>-a);y=y(y<a);   % Truncate
[~,f_Levy,yi,~]=kde(y(1:1:end)',2^16,min(y)-0.1*range(y),max(y)+0.1*range(y));
plot(yi,f_Levy,'-b');xlim([-3 3]);ylim([0 0.5]);hold on;

L1=-4;R1=4;
mesh=linspace(L1,R1,2000);
D1=-mesh.^3;
D2=0.*mesh+0.5;DD2=0.*mesh;
tf=50000;   % tf is the final time we use to solve Fokker-Plank equation to get stationar density (in theory it should be infinite) 
si=0.1;
PDF0=@(x)1/sqrt(2*pi*si^2)*exp(-0.5*((x-mean([L1 R1]))/si)^2); % An initial Gaussian distribution centered at the state space and with std of si=0.1
f_Levy=FP(L1,R1,mesh,D1,D2,DD2,PDF0,tf);  % FP numerically solves the Fokker-Planck equation (only for Gaussian white noise)
plot(mesh,f_Levy,'-k');xlim([-3 3]);ylim([0 0.5]);
xlabel('System state');ylabel('Stationary distribution');hold on;title('Quartic potential (monostable)');

%%
% Figure 3 (Reconstruction of two datasets with the same density)
clc;

% You can generate new data as bellow
% dt=10^(-2);T=10^5;x=zeros(1,T);y=x;
% for i=2:T
%     x(i)=x(i-1)-x(i-1)*dt+sqrt(2*dt)*randn;
%     y(i)=y(i-1)+(y(i-1)-y(i-1)^3)*dt+sqrt(2*dt*(1+y(i-1)^2))*randn;
% end
% data1=x(1:end);data2=y(1:end);

% A (The data and their distribution)
data1=readmatrix('FirstDataset.csv');data2=readmatrix('SecondDataset.csv');
l=linspace(-4,4,2000);h=[0 0];
[f1,~]=ksdensity(data1,l);f=1./sqrt(2.*pi).*exp(-0.5.*l.^2);
figure,plot(data1,'-k');ylim([-4 4]);set(gca,'xtick',[]);set(gca,'ytick',[]);
figure,h(2)=plot(l,f1,'-r','LineWidth',3);hold on;h(1)=plot(l,f,'-.k','LineWidth',3);set(gca,'xtick',[]);set(gca,'ytick',[]);ylim([0 0.43]);
set(gca,'xtick',[]);set(gca,'ytick',[]);

l=linspace(-4,4,2000);
[f2,~]=ksdensity(data2,l);f=1./sqrt(2.*pi).*exp(-0.5.*l.^2);
figure,plot(data2,'-k');ylim([-4 4]);set(gca,'xtick',[]);set(gca,'ytick',[]);
figure,h(2)=plot(l,f2,'-r','LineWidth',3);hold on;h(1)=plot(l,f,'-.k','LineWidth',3);ylim([0 0.43]);
set(gca,'xtick',[]);set(gca,'ytick',[]);

% B&C&D&E
Tau=1:5;method='Nadaraya-Watson1';L=-3;R=3;bins=50;dt=10^(-2);
[D1,D2,D4,ErrorD1,ErrorD2,ErrorD4,N,C]=LangevinReconst(data1,L,R,bins,Tau,dt,method);
[D11,D22,D44,ErrorD11,ErrorD22,ErrorD44,N,C]=LangevinReconst(data2,L,R,bins,Tau,dt,method);

h=[0 0];
figure,
h(1)=plot(C,-C,'--r','LineWidth',1.5);hold on;h(2)=plot(C,D1,'.k');hold on;errorbar(C,D1,ErrorD1,'.k');hold on;plot(C,0.*C,'-k');
xlabel('System state');ylabel('Rate of change');ylim([-3 3]);
xticks([-3 0 3]);xticklabels([-3 0 3]);yticks([-3 0 3]);yticklabels([-3 0 3]);
legend(h,{'True','Estimated'});
figure,
h(1)=plot(C,0.*C+1,'--r','LineWidth',1.5);hold on;h(2)=plot(C,D2,'.k');hold on;errorbar(C,D2,ErrorD2,'.k');ylim([0 inf]);
xlabel('System state [a.u.]');ylabel('Noise intensity');
xticks([-3 0 3]);xticklabels([-3 0 3]);yticks([0 0.6 1.2]);yticklabels([0 0.6 1.2]);
legend(h,{'True','Estimated'});
figure,
h(1)=plot(C,C-C.^3,'--r','LineWidth',1.5);hold on;h(2)=plot(C,D11,'.k');hold on;errorbar(C,D11,ErrorD11,'.k');hold on;plot(C,0.*C,'-k');
yticks([-25 0 25]);yticklabels([-25 0 25]);xticks([-3 0 3]);xticklabels([-3 0 3]);
xlabel('System state [a.u.]');ylabel('Rate of change');ylim([-25 25]);
legend(h,{'True','Estimated'});
figure,
h(1)=plot(C,1+C.^2,'--r','LineWidth',1.5);hold on;h(2)=plot(C,D22,'.k');hold on;errorbar(C,D22,ErrorD22,'.k');
xlabel('System state [a.u.]');ylabel('Noise intensity');ylim([0 inf]);
yticks([0 5 9]);yticklabels([0 5 9]);xticks([-3 0 3]);xticklabels([-3 0 3]);
legend(h,{'True','Estimated'});

%%
% Figure 4 (The Leavy reconstruction of synthetic data)

% Note: This section is going to take around 3 min
clc;
S=load('LandauData.mat');data=S.data;
L=-2.;R=2.2;  % Here, we have considered a smaller data range to make sure all bins receive enough data (this, usually, needs a bit of experimentation).
% data(data<L)=nan;data(data>R)=nan;
Tau=1:5;bins=50;dt=10^(-2);
method='Nadaraya-Watson';
Iter=4;  % this reconstruction algorithm is iterative. Normally 2 iterations are enough but we use 4
useparfor=1;  % set this to 0 if you do not have parallell toolbox
[g,Alpha,h,Errorg,stdAlpha,Errorh,N,C]=LevyReconst(data,L,R,bins,Tau,dt,method,Iter,useparfor);

figure,plot(data,'-k');ylim([L R]);yticks([-2 0 2]);yticklabels([-2 0 2]);set(gca,'xtick',[]);
hh=[0 0];
figure,
hh(1)=plot(C,C-C.^3,'--r','LineWidth',1.5);hold on;hh(2)=plot(C,g,'.k');hold on;errorbar(C,g,Errorg,'.k');hold on;plot(C,0.*C,'-k');
xlim([L R]);ylim([-inf inf]);xticks([-2 -1 0 1 2]);yticklabels([-2 -1 0 1 2]);yticks([-6 -3 0 3 6]);yticklabels([-6 -3 0 3 6]);
legend(hh,{'True','Estimated'});
figure,
hh(1)=plot(C,0.*C+0.3,'--r','LineWidth',1.5);hold on;hh(2)=plot(C,h,'.k');hold on;errorbar(C,h,Errorh,'.k');
xlim([L R]);ylim([0 inf]);xticks([-2 -1 0 1 2]);yticklabels([-2 -1 0 1 2]);yticks([0 0.3 0.6]);yticklabels([0 0.3 0.6]);ylim([0 0.6]);
legend(hh,{'True','Estimated'});

%%

% Figure 5 (The Leavy reconstruction of climate data)

% Note: This section is going to take around 3 min
clc;
data=csvread('Cadata.csv');data = fliplr(data);L=-2+0.3;R=max(data)-0.2;%data=data(data>=L & data<=R);
data(isinf(data))=nan;   % note that this dataset has 43 inf values (they must turn into nan).

Tau=1:5;bins=50;dt=10^(-2);
method='Nadaraya-Watson';
Iter=4;  % this reconstruction algorithm is iterative. Normally 2 iterations are enough but we use 4
useparfor=1;  % set this to 0 if you do not have parallell toolbox
[g,Alpha,h,Errorg,stdAlpha,Errorh,N,C]=LevyReconst(data,L,R,bins,Tau,dt,method,Iter,useparfor);

w_g=1./Errorg;w_h=1./Errorh;
[fit_g, ~] = FitSmoothinSpline(C, g, w_g);
[fit_h, ~] = FitSmoothinSpline(C, h, w_h);

figure,plot(data,'-k');ylim([L R]);yticks([-1.5 0 1.5 3]);yticklabels([-1.5 0 1.5 3]);set(gca,'xtick',[]);
hh=[0 0];
figure,
hh(1)=plot(C,g,'.k');hold on;errorbar(C,g,Errorg,'.k');hold on;plot(C,0.*C,'-k');hold on;hh(2)=plot(C,fit_g(C),'--r','LineWidth',1.5);
xlim([L R]);ylim([-inf inf]);xticks([-1.5 0 1.5 3]);xticklabels([-1.5 0 1.5 3]);yticks([-12 -6 0 6]);yticklabels([-12 -6 0 6]);
legend(hh,{'Estimated','Fited curve'});
figure,
hh(1)=plot(C,h,'.k');hold on;errorbar(C,h,Errorh,'.k');hold on;hh(2)=plot(C,fit_h(C),'--r','LineWidth',1.5);
xlim([L R]);xticks([-2 -1 0 1 2]);xticks([-1.5 0 1.5 3]);xticklabels([-1.5 0 1.5 3]);yticks([0 0.5 1 1.5]);yticklabels([0 0.5 1 1.5]);ylim([0 1.5]);
legend(hh,{'Estimated','Fited curve'});

%%
% Figure 6 (Exit time)
clc;
syms x mu1(x) sigma1(a) mu2(x) sigma2(a) 
mu1(x)=-x;sigma1(x)=sqrt(2);mu2(x)=(x-x^3);sigma2(x)=sqrt(2*(1+x^2));
domain=[-2 2];BC='AA';

T1=MeanExitTime('Parametric',mu1,sigma1,domain,BC);
T2=MeanExitTime('Parametric',mu2,sigma2,domain,BC);

h=[0 0];
h(1)=plot(T1,'-k');hold on;h(2)=plot(3.*T2,'-r');ylim([-inf 5.5]);
legend(h,{sprintf('First datase in Figure 3'), sprintf('Second datase in Figure 3')});
xticks([-2 -1 0 1 2]);xticklabels([-2 -1 0 1 2]);yticks([0 1 2 3 4 5]);yticklabels([0 1 2 3 4 5]);
title('Exit time as a measure of resilience');xlabel('System state');ylabel('Mean exit time');

%%
% Figure 7(Tree cover data)

clc;
% Calculating the density, density derivative, the dynamivs, and estimating the noise intensity

load('DatSA.mat');
R=100;
p=1650;

I1=find(DatSA(:,1)==p-25,1);I2=find(DatSA(:,1)==p+26,1);I2=I2-1;
data=DatSA(I1:I2,2:end);
data=data+(0.5.*rand(size(data))).*(data<=0.5)+(rand(size(data))-0.5).*(data>0.5)+(-0.5.*rand(size(data))).*(data>=99.5); % This is to perturb the data a bit to smooth them a bit

data1=reshape(data,[1,size(data,1)*size(data,2)]);
data2=[data nan(size(data,1),1)];   % I put a nan at the end of each data at all spatial locations
data2=reshape(data2,[1,size(data2,1)*size(data2,2)]);
X0=data2(1:end-1);X=data2(2:end);

if length(data1)>3*10^6
    data1=randsample(data1,3*10^6);
end
mesh=linspace(0,max(data1),3000);
[f,~]=ksdensity(data1,mesh);
F=@(x)interp1(mesh,f,x,'cubic');
pp_f=interp1(mesh,f,'cubic','pp');
pp_df=fnder(pp_f);
dF=@(x)ppval(pp_df,x);
dyn=@(x)interp1(mesh,dF(mesh)./F(mesh),x,'cubic');
dt=1;  % this is arbitrary (although the outcome of reconstruction depends on dt. For instance, if you choose dt=0.01 then the final result will be 100 times bigger)
cost=@(sigma)Cost(sigma,dyn,X0,X,dt);

%For simulations
% T=10^4;x=zeros(1,T);SIGMA=12.4;
% for i=2:T
%     x(i)=x(i-1)+(SIGMA^2/2*dyn(x(i-1)))*dt+sqrt(dt)*SIGMA*randn;
%     if x(i)<=0
%         x(i)=dt;
%     end
% end

lb=0;    % lower bound for the single unknown parameter (noise intensity). 
ub=1000; % upper bound (based on the chosen 'dt', the upper bound should be big enough)
SIGMA=fminbnd(cost,lb,ub);   % SIGMA is our final goal, i.e., noie intensity


% Plot the spatial locations

% addpath('C:\Users\P264168\Desktop\TreeCover(for ASS)');
% p1=1600;p2=1650;
% load('map_trp.mat');load('mod00trp.mat');load('human_trp.mat');pr=map_trp;
% T=20000; prSA=pr(:,1:T);
% ISA=double(h5read('ISA.h5', '/DS1'));
% I=prSA>=p1 & prSA<=p2;
% % Itreeless=ISA & prSA>=0;
% Iocean=prSA<0;
REGIONS=zeros(T,T);REGIONS(Iocean==0)=5;REGIONS(I==1)=3;imagesc(REGIONS);
set(gca,'xtick',[]);set(gca,'ytick',[]);

h=[0 0];
h(1)=plot(mesh,SIGMA.^2./2.*dyn(mesh),'-k','LineWidth',1);hold on;h(2)=plot(mesh,1150.*f,'-.b');hold on;plot(mesh,0.*mesh,'-k');
ylim([-20 40]);xlim([-inf inf]);legend(h,{'Determinist forces','Distribution'});
xticks([0 30 60 90]);xticklabels([0 30 60 90]);yticks([-20 0 20 40]);yticklabels([-20 0 20 40]);

%**************************
%******Subroutine**********
%**************************
function cost=Cost(sigma,dyn,X,X0,dt)
mu_X0=sigma^2/2*dyn(X0);
LL=-1./2.*(log(2.*pi.*dt.*sigma.^2)+((X-X0-mu_X0.*dt)./(sigma.*sqrt(dt))).^2);
cost=double(-sum(LL,'omitnan'));   %objective 
end









    


