function dwtr = dwtr(data, L, filterh)
%  function dwtr = dwtr(data, L, filterh); 
%  Calculates the DWT of periodic data set
%  with scaling filter  filterh  and  L  detail levels. 
%
%   Example of Use:
%   data = [1 0 -3 2 1 0 1 2]; filterh = [sqrt(2)/2 sqrt(2)/2];
%   wt = dwtr(data, 3, filterh)
%--------------------------------------------------------------------------

n = length(filterh);                
C = data(:)';                             
dwtr = [];                                                           
H = filterh;                          
G = fliplr(filterh);
G(1:2:n) = -G(1:2:n);           
for j = 1:L                              
   nn = length(C);                  
   C = [C(mod((-(n-1):-1),nn)+1)  C]; 
   D=filter(G,[1],C);                                  
   D = D([n:2:(n+nn-2)]+1);     
   C=filter(H,[1],C);                                     
   C = C([n:2:(n+nn-2)]+1);     
   dwtr = [D,dwtr];                                           
end;                                    
dwtr = [C, dwtr];               

