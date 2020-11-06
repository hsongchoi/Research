function  data = idwtr(wtr, L, filterh)
% function data = idwt(wtr, L, filterh); Calculates the IDWT of wavelet
% transformation wtr using wavelet filter  "filterh"  and  L  scales.  
% Use
%>>  t=linspace(0,1,1024); data = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
%>> filter=[sqrt(2)/2 sqrt(2)/2];
%>> max(abs(data - idwtr(dwtr(data,3,filter), 3,filter)))
%
%ans = 3.8858e-016

nn = length(wtr);   n = length(filterh);  
if nargin==2, L = round(log2(nn)); end;        
  H = fliplr(filterh);                                   
  G = filterh; G(2:2:n) = -G(2:2:n);
LL = nn/(2^L);                                   
C =  wtr(1:LL);                                
for j = 1:L                                      
   w  = mod(0:n/2-1,LL)+1;                      
   D  = wtr(LL+1:2*LL);                       
   Cu(1:2:2*LL+n) = [C C(1,w)];                   
   Du(1:2:2*LL+n) = [D D(1,w)];                    
   C  = filter(H,[1],Cu) + filter(G,[1],Du); 
   C  = C([n:n+2*LL-1]-1);                       
   LL = 2*LL;                                      
end;
data = C;                                         

