
function [g,Alpha,h,Errorg,stdAlpha,Errorh,N,C]=LevyReconst(data,L,R,bins,Tau,dt,method,Iter,useparfor)

% This code reconstructs a one-dimensional Levy model dx = g(x)dt+h(x) dL based on the following paper

% Siegert, S., & Friedrich, R. (2001). Modeling of nonlinear Lévy processes by data analysis. Physical Review E, 64(4), 041107.

% where dL reperesents an alpha-stable levy motion and in general has 4 pasrameters: 1) alpha (stability index, a number in the interval
% [0 2]) 2) betta (skewness, a number in [-1 1]) 3)  dispersion or scale parameter 4)center. Here, dispersion parameter is set to 1 (since
% the function h(x) can take responsibility for it), center is set to 0 (since the function h(x) can take responsibility for it), skewness
% parameter is 0 since the reference paper only considers this case (i.e., symmetric Levy motion). The goal of this algorithm, therefore, is
% to estimate the functions g(x) and h(x), and the stability index.
% If stability index is 2 then the Levy motion dL reduces to the standard Brownian motion dW. The smaller the stability index is the jumpier the data would be.

% Note 1: the methodology in the above reference works when stability index is in the interval (1,2].
% Note2: This code can only work when stability index is in the interval [1.2, 2]. The reason has to do with our limited capacity in estimating
% some pre-calculated quantities (otherwise, this code would be extremely slow) which are computationaly very expensive. Estimation of those pre-calculates
% quntities become computationally intensive as alpha aproaches 1. Nevertheless, in most applications this range is enough.
% Note 3: the dataset can contain NaN values and this code can handle it.

%Inputs:
    %   data       the data set (evenly sampled time series) 
    %   L and R    are, respectively, the lower and upper limit considered for the states being spanned by the range
    %              of data. Normally, L=min(data) and R=max(data) are good choices but for limited data one might consider a slightly
    %              bigger value for L and a slightly smaller value for R (to make sure all bins receive enough, say 100, number of observations).
    %   bins       number of bins for the reconstruction Should be carefully chosen depending on the size of data. This number should not be so big so that
    %              there would be little data in each bin and it should not be so small so that we are left with few bins. It is recommended to choose
    %              bins in such a way that each bin will, at least, contain 100 observations).
    %   Tau        a list of time lags considered (integer values) (we chose 1:5. For high resolution data we can consider bigger number of time lags).
    %   dt         time step of the data series (for real datasets feel free to choose any positive number you like)
    %   method     is either 'Nadaraya-Watson' or empty, 'Nadaraya-Watson' is recommended
    %              as the result is smoothed.
    %   Iter       this reconstruction algorithm is an iterative algorithm. In the above citation it is mentioned that 2 iterations is enough 
    %              (This is tested by I. It is recommended to choose 3-4 iterations)
    %   useparfor  set it to 1 if you have parallel toolbox, otherwise set it to 0
    %
%Outputs:
    %   g         vector of estimated deterministic function for each bin 
    %   Alpha     estimated stability index
    %   h         vector of estimated stochastic function (noise intensity) for each bin 
    %   Errorg    vector of estimated error in g expressed as standard deviation for each bin.
    %   stdAlpha  standard deviation of estimated stability index Alpha
    %   Errorh    vector of estimated error in h expressed as standard deviation for each bin.
    %   N         vector representing the number of data per bin.
    %   C         vector of bin centers (So, if you want to plot g or h functions then the the proper commands
    %             are plot(C,g,'.k') and plot(C,h,'.k')), respectively.
    %
    %   this function requires the Curve Fit Toolbox to be installed
    %   this function requires the Chebfun package. You can find the latest version of Chebfun in the following link:
    %   https://nl.mathworks.com/matlabcentral/fileexchange/47023-chebfun-current-version
    %
    %   Implemented in Matlab by Babak M.S. Ariani

if useparfor
    numWorkers = Inf;
else
    numWorkers = 1;
end

% Here, we benifit from some pre-calculated quntities (otherwise, the code could take hours)
F_alpha=readmatrix('Fa3(20000Samples).csv');F_alpha=F_alpha(21:end);F1=F_alpha;
l=1.2:0.01:2;StartPoint=ones(1,11);[fit,~]=createFitRational(l,F1,'rat55',StartPoint);
Fa=chebfun(fit(linspace(1.2,2,2000)),[1.2 2],'equi',1000);dFa=diff(Fa);
Delta70=readmatrix('Delta70.csv');Delta70=Delta70(21:end);Delta70=chebfun(Delta70',[1.2 2],'equi',2000);

TAU=dt.*Tau;

% In bellow, I bin the data based on the pre-specified lower and upper limits L and R. However, in the process of 
% reconstruction all the data are participated
data1=data;
x=linspace(L,R,bins+1);                % x is the bins borders.
dx=x(2)-x(1);C=x+dx/2;C=C(1:end-1);    % C is the bin centers.
[N,~]=histcounts(data1,x);             % N is the amount of data in each bin.

% This is an interative algorithm. In the first iteration we assume unit uncerainty for all unknowns. This will be updated in the next iterations
uncertaintygY_axis=ones(bins,length(Tau));uncertaintyAlphaY_axis=ones(bins,length(Tau));uncertaintyhY_axis=ones(bins,length(Tau));

g=zeros(1,bins);Errorg=zeros(1,bins);   % g is the deterministic force.
alpha=zeros(1,bins);                    % alpha is the stability index
h=zeros(1,bins);Errorh=zeros(1,bins);   % h is the stochastic force (noise intensity).

A1=zeros(1,bins);T1=zeros(length(Tau),bins);
A2=zeros(1,bins);T2=zeros(length(Tau),bins);

for ITER=1:Iter
    if strcmp(method,'Nadaraya-Watson')==0
        for tau=Tau   % This is the direct way to estimate conditional methods.
            for i=1:bins
                I=find(data>=x(i) & data<x(i+1));
                I=setdiff(I,I(I+tau>length(data))); % Clearly, data shifted in time by tau which are not in the time domain must be excluded from the analysis.
                Iplus=I+tau;                       % We are now sure that Iplus does not exceed the maximum time domain of data.
                A1(i)=mean(data(Iplus)-data(I),"omitnan");
            end
            T1(tau,:)=A1;
        end
    else
        % Nadarya-Watson estimation of conditional moments
        % hh=(4/(3*length(dat)))^(1/5)*std(dat);K=@(x)exp(-x.^2./2);   %Gaussian kernel (If you wish to use this)
        hh=1.048*(length(data))^(-1/5)*std(data,"omitnan");K=@(x)(x>-sqrt(5) & x<sqrt(5)).*3.*sqrt(5)./100.*(5-x.^2);  %Epanechnikov kernel
        for tau=Tau
            parfor(i = 1:bins, numWorkers)
            % for i=1:bins
                DATA=data;   % dat is a broadcast variablle
                data1=DATA(1:tau:end);plot(data1,'-k');
                A=K((data1(1:end-1)-C(i))./hh);B=diff(data1);
                N1=sum(A.*B,"omitnan");E1=sum(A,"omitnan");  % Note that I use N1 to distinguish it with N above.
                A1(i)=N1/E1;
            end
            T1(tau,:)=A1;
        end
    end

    % Estimating the deterministic force (the function g(x))
    for i=1:bins
        [fit,~]=LinearFit(TAU,T1(:,i)',uncertaintygY_axis(i,:));Errorg(i)=SlopeError(TAU,T1(:,i)',fit(TAU)',uncertaintygY_axis(i,:));
        g(i)=differentiate(fit,0);
    end

    if strcmp(method,'Nadaraya-Watson')==0
        for tau=Tau   % This is the direct way to estimate conditional methods.
            for i=1:bins
                I=find(data>=x(i) & data<x(i+1));
                I=setdiff(I,I(I+tau>length(data)));
                Iplus=I+tau;
                A2(i)=mean(abs(data(Iplus)-data(I)-g(i)*TAU(tau)),"omitnan");
            end
            T2(tau,:)=A2;
        end
    else
        % Nadarya-Watson estimation of conditional moments
        for tau=Tau
            parfor(i = 1:bins, numWorkers)
                DATA=data;TAU1=TAU;   % dat and TAU are broadcast variables
                data1=DATA(1:tau:end);
                A=K((data1(1:end-1)-C(i))./hh);B=abs(diff(data1)-g(i)*TAU1(tau));
                N1=sum(A.*B,"omitnan");E1=sum(A,"omitnan");
                A2(i)=N1/E1;
            end
            T2(tau,:)=A2;
        end
    end
    
    % Estimating the stability index (alpha)
    % a=2.^(0:length(Tau));talpha=a(a<=max(Tau));Talpha=dt.*talpha;   %
    % Note that to estimate alpha we could use time steps of the form dt,2dt,2^2dt,...2^n*dt as is mentioned in the reference paper but i is not necessary (if you wish to do so, uncomment this line and the first line after the next for-loop, but comment the line after that)
    parfor(i = 1:bins, numWorkers)
        %   [fit,gof]=createFit(log(Talpha),log(T2(talpha,i))',ones(1,length(talpha)),fitform,opts);   %Based on the above reference where 2^n time steps are used. I did not follow this (if you wish to do so, then uncomment this line and comment the next line).
        [fit,~]=LinearFit(log(TAU),log(T2(:,i))',uncertaintyAlphaY_axis(i,:));
        alpha(i)=1/differentiate(fit,0);    % alpha is the inverse of the fited line to ln(T2) over ln(tau).
    end
    
    idx=alpha<2 & alpha>1.2;
    Alpha=[];stdAlpha=[];
    if ~isempty(idx)
        Alpha=mean(alpha(idx));stdAlpha=std(alpha(idx));
    else
        disp('This dataset has a very bad quality so that stability index for all bins falls outside the legitimate range of [1 2]');
        return;
    end

    % Estimating the stochastic force (the function h(x))
    parfor(i = 1:bins, numWorkers)
        Fa1=Fa;   % Fa is a broadcast variable
        [fit,~]=LinearFit(TAU,T2(:,i)'./(TAU.^(1./Alpha-1).*Fa1(Alpha)),uncertaintyhY_axis(i,:));
        Errorh(i)=SlopeError(TAU,T2(:,i)'./(TAU.^(1./Alpha-1).*Fa1(Alpha)),fit(TAU)',uncertaintyhY_axis(i,:));
        h(i)=differentiate(fit,0);
    end

    %*****************************************************
    %Uncertainty analysis (updating the old uncertainties)
    %*****************************************************

    %All uncertainties are state-by-lags matrix.

    %Uncertainty of g (This is, indeed, the uncertainty of the numerator (T1) on the Y_axis).
    n=length(TAU);
    uncertaintygY_axis=Delta70(Alpha).*repmat(h',[1 n])./repmat(N',[1 n]).^(1-1./Alpha).*repmat(TAU,[bins 1]).^(1./Alpha);

    %Uncertainty of h (This is, indeed, uncertainty of tau*T3).
    delABS=zeros(1,bins);
    parfor(i = 1:bins, numWorkers)
        delABS(i)=Delta70ABS(Alpha,N(i),Fa);
    end
    delABS=2.*delABS;   % I did not understand why in my calculations delABS is half those in Silke paper. So, this is whay I multiply it by 2! But, I noticed that this has very little effect on the final results and hard to see by eye.
    A=repmat(h',[1 n]).*repmat(delABS',[1 n])./Fa(Alpha).*repmat(TAU,[bins 1]);
    B=T2'.*repmat(TAU,[bins 1]).^(1-1./Alpha).*1./Fa(Alpha).*abs(dFa(Alpha)./Fa(Alpha)-repmat(log(TAU),[bins 1])./Alpha.^2).*stdAlpha;
    uncertaintyhY_axis=A+B;

    %Uncertainty of 1/Alpha
    uncertaintyAlphaY_axis=repmat(h',[1 n])./T2'.*repmat(delABS',[1 n]).*repmat(TAU,[bins 1]).^(1./Alpha);
end
idx=find(alpha>2);
if ~isempty(idx)
disp(['At bins {' num2str(idx) '} stability index bigger than 2 is obtained. If such bins are close to data borders then this is not a big issue, otherwise consider shrinking the data domain. Note also that the code ignores such bins when it calculates the stability index']); 
end

end


