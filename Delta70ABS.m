function D=Delta70ABS(alpha,N,Fa)

k=5000;R=10;
Z=mean(reshape(abs(alpha_stable_sym(alpha,k*N)),[N k]));Z=Z(Z<=R);
[f,xi]=ksdensity(Z,linspace(min(Z),R,1000));f=chebfun(f',[min(Z) R],'equi',1000);
Cost=@(u)abs(sum(f,[Fa(alpha)-u Fa(alpha)+u])-0.7);
lb=0;ub=min(Fa(alpha)-min(Z),R-Fa(alpha));x0=1;

% fminbnd optimization
options=optimset('TolX',1e-6);
D=fminbnd(Cost,lb,ub,options);

% pattersearch optimization
% options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'UseVectorized',false,'TolX',1e-6);
% D=patternsearch(Cost,x0,[],[],[],[],lb,ub,[],options);


end

