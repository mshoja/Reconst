
function [D1,D2,D4,ErrorD1,ErrorD2,ErrorD4,N,C]=LangevinReconst1(data,L,R,bins,Tau,dt,method)

 %Inputs:
    %   data      the data set (evenly sampled time series) (can also take
    %             segmented data (see makesegments)
    %   L and R   are, respectively, the lower and upper limit considered for the states being spanned by the range
    %             of data. Normally, L=min(data) and R=max(data) are good choices but for limited data one might consider a slightly
    %             bigger value for L and a slightly smaller value for R.
    %   bins      number of bins for the reconstruction Should be carefully chosen depending on the size of data. This number should not be so big so that
    %             there would be little data in each bin and it should not be so small so that we are left with few bins. In the following citation, it is recommended to choose
    %             bins in such a way that each bin will, at least, contain 100 observations).
    %   Tau       a list of time lags considered (integer values) (we chose 1:5. For high resolution data we can consider bigger number of time lags).
    %   dt        time step of the data series
    %   method    is either 'Nadaraya-Watson' or empty, 'Nadaraya-Watson' is recommended
    %             as the result is smoothed.
    %
    %Outputs:
    %   D1        vector of estimated drift function for each bin (first Kramers-Moyal coefficient)
    %   D2        vector of estimated diffusion function for each bin (second Kramers-Moyal coefficient)
    %   D4        vector of estimated fourth Kramers-Moyal coefficient for each bin (should be small)
    %   ErrorD1   vector of estimated error in D1 expressed as standard deviation for each bin.
    %   ErrorD2   vector of estimated error in D2 expressed as standard deviation for each bin.
    %   N         vector representing the number of data per bin.
    %   C         vector of bin centers (So, if you want to plot drift or diffusion functions then the the proper commands
    %             are plot(C,D1,'.k') and plot(C,D2,'.k')), respectively.
    %
    %   this function requires the Curve Fit Toolbox to be installed
    %
    %   Implemented in Matlab by Babak M.S. Ariani
    %   You can alternatively use the R package with instructions about the method in
    %     Rinn, P., Lind, P. G., Wächter, M. & Peinke, J. The Langevin Approach: An R Package for Modeling Markov Processes. arXiv preprint arXiv:1603.02036 (2016).


% In bellow, I bin the data based on the pre-specified lower and upper limits L and R. However, in the process of 
% reconstruction all the data are participated
data1=data;
x=linspace(L,R,bins+1);                % x is the bins borders.
dx=x(2)-x(1);C=x+dx/2;C=C(1:end-1);    % C is the bin centers.
[N,~]=histcounts(data1,x);             % N is the amount of data in each bin.

% Estimating conditional moments
M1=zeros(1,bins);M2=M1;M4=M1;   
A1=zeros(length(Tau),bins);A2=A1;A4=A1;    % Ai Matrix contain the Mi vectors for different values of Tau. 
ErrorM1=zeros(length(Tau),bins);ErrorM2=ErrorM1;     % Error 1 and Error2 are variance of estimmated M1 and M2 coefficients for different lags and bins.

% This is the direct way to estimate conditional methods.
if strcmp(method,'Nadaraya-Watson')==0  
    for tau=Tau
        for i=1:bins
            I=find(data>=x(i) & data<x(i+1));I=setdiff(I,I(I+tau>length(data)));Iplus=I+tau;
            M1(i)=mean(data(Iplus)-data(I),"omitnan");M2(i)=mean((data(Iplus)-data(I)).^2,"omitnan");M4(i)=mean((data(Iplus)-data(I)).^4,"omitnan");
        end
        A1(tau,:)=M1;A2(tau,:)=M2-M1.^2;A4(tau,:)=M4;  % Note, to account for finite-tau corrections to diffusion consider A2(tau,:)=M2-M1.^2.
        ErrorM1(tau,:)=(M2-M1.^2)./N;ErrorM2(tau,:)=(M4-M2.^2)./N;
    end
else
% Nadarya-Watson estimation of conditional moments
h=1.048*(length(dat))^(-1/5)*std(dat);K=@(x)(x>-sqrt(5) & x<sqrt(5)).*3.*sqrt(5)./100.*(5-x.^2);  %Epanechnikov kernel
% h=(4/(3*length(dat)))^(1/5)*std(dat);K=@(x)exp(-x.^2./2); % In case you wish to use Gaussian Kernel.
for tau=Tau 
    parfor i=1:bins    % Here, I use Matlan parallel computing toolbox. If you do not have it just use 'for loop' not 'parfor loop'.
        data1=data(1:tau:end);
         A=K((data1(1:end-1)-C(i))./h);B=diff(data1);
         N1=sum(A.*B,"omitnan");N2=sum(A.*B.^2,"omitnan");N4=sum(A.*B.^4,"omitnan");E=sum(A,"omitnan");
         M1(i)=N1/E;M2(i)=N2/E;M4(i)=N4/E;
    end
    A1(tau,:)=M1;A2(tau,:)=M2-M1.^2;A4(tau,:)=M4;
    ErrorM1(tau,:)=(M2-M1.^2)./N;ErrorM2(tau,:)=(M4-M2.^2)./N;
end
end

WeightM1=1./ErrorM1;
WeightM2=1./ErrorM2;

% This is to make sure all weights are positive.
for i=1:bins
    if nnz(WeightM1(:,i)<0)>0
        WeightM1(:,i)=ones(1,length(Tau));
    end
end
for i=1:bins
    if nnz(WeightM2(:,i)<0)>0
        WeightM2(:,i)=ones(1,length(Tau));
    end
end

% Estimating the limit conditional moments/dt as dt-->0
D1=zeros(1,bins);D2=D1;D4=D1;Tau=dt.*Tau;
ErrorD1=D1;ErrorD2=D1;ErrorD4=D1;

% Calculation of Drift (D1), Diffusion (D2) and, D4 
for i=1:bins
    
    [fit1,gof1]=LinearFit(Tau,A1(:,i)',WeightM1(:,i)');ErrorD1(i)=SlopeError(Tau,A1(:,i)',fit1(Tau)',WeightM1(:,i)');
    [fit2,gof2]=LinearFit(Tau,A2(:,i)',WeightM2(:,i)');ErrorD2(i)=SlopeError(Tau,A2(:,i)',fit2(Tau)',WeightM2(:,i)');  
    [fit4,gof4]=LinearFit(Tau,A4(:,i)',ones(1,length(Tau)));ErrorD4(i)=SlopeError(Tau,A4(:,i)',fit4(Tau)',ones(1,length(Tau)));
    
    D1(i)=differentiate(fit1,0);   
    D2(i)=1/2*differentiate(fit2,0);
    D4(i)=1/24*differentiate(fit4,0); 
end
end
