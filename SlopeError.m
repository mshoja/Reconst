function e=SlopeError(x,y,y_est,w)
n=length(x);
% w=w./sum(w);
e=sqrt(sum(w.*(y-y_est).^2)/((n)*(sum(w.*(x-mean(x)).^2))));
end