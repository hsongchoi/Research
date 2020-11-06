function output=ker_sm(data, sigma);

N=length(data);
[B,A]=ndgrid(data,data);
output=sum(normpdf(B-A,0,sigma))/N;
