function output=b_factor(d,mu,tau);
num= ((2*mu)^0.5)/2 .* exp(-(2*mu)^0.5 .* abs(d)) .* (2*tau^2-1/mu) ;
denom=(tau*exp(-abs(d)./tau) -exp(-(2*mu)^0.5 .* abs(d))./sqrt(2 * mu));
output=num./denom;
