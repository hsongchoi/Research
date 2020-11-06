function output=b_factor(d,h,m)
if (-m > d)
    num= h/2 .* exp(-h .* abs(d));
denom=(exp(h.*d).*sinh(h.*m))./(2*m) ;
 output=num./denom;
elseif  (d>= -m) & (m >= d)
       num= h/2 .* exp(-h .* abs(d));
denom=1/(2*m)-exp(-h.*m).*cosh(h.*d)./(2*m);
 output=num./denom;
else 
     num= h/2 .* exp(-h .* abs(d));
denom=(exp(-h.*d).*sinh(h.*m))./(2*m);
 output=num./denom;
end
    

