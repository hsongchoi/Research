function dwtr = dwtrn(data, L, filterh)
%  function dwtr = dwtr_n(data, L, filterh); 
%  Calculates the DWT of periodic data set (multi-dimensional)
%  with scaling filter  filterh and L detail levels.
%--------------------------------------------------------------------------

dwtr = data;
d = ndims(dwtr);                                   % number of dimensions
ssize = size(dwtr);                                % size of the space

oned = ( d==2 && ~isempty(find(ssize==1)) );       % true, if vector (row or column)
    
if(oned)  msize = max(ssize);
else      msize= min(ssize);  end;

L = min([L floor(log2(msize))]);
for i = 1:L,                                       % for each level       
    sindex = get_space_index(ssize);               % subspace index 
    eval(['dwtr(' sindex ') = dwtr_n1(dwtr(' sindex '), filterh);']);    % replace subspace
    ssize = floor(ssize/2);                        % new subspace size
    if(oned) ssize(find(ssize==0)) = 1; end
end
                                  

function dwtr = dwtr_n1(data, filterh)
%  function dwtr = dwt_n1(data, filterh); 
%  Calculates the DWT of periodic data set (n dimensional)
%  with scaling filter  filterh and 1 detail level 
%--------------------------------------------------------------------------

d = ndims(data);                                   % dimension of the data
s = size(data);                                    % size of each dimension

oned = ( d==2 && ~isempty(find(s==1)) );           % true, if vector (row or column)

if ( oned )                                        % one dimensional data, end of recursion            
    C_r = reshape(data, [length(data) 1]);         % reshape to one dimensional column vector
    C_t = dwtr_11(C_r, filterh);                   % call one dimensional specific function    
    dwtr = reshape(C_t, size(data));               % reshape back to original shape
    return;
end

dwtr = data;                             
for i=1:d,                                                    % for each dimension
    for j=1:s(i),                                             % for each slice within dimension 
        C_s = eval(['dwtr(' get_slice_index(d, i, j) ')']);   % get the slice
        if (d==2)                                             % if the slice is one dimensional 
            C_r = reshape(C_s, [length(C_s) 1]);              % reshape to column vector    
        else                                                   
            C_r = squeeze(C_s);                               % remove the virtual dimension
        end
        C_t = dwtr_n1(C_r , filterh);                         % solve the n-1 dimensional problem 
        C_r = reshape(C_t, size(C_s));                        % reshape back to slice dimensions
        eval(['dwtr(' get_slice_index(d, i, j) ') = C_r;']);  % set old slice to new slice  
    end
end


function dwtr = dwtr_11(data, filterh)
%  function dwtr = dwt_11(data, filterh); 
%  Calculates the DWT of periodic data set (row vector)
%  with scaling filter  filterh  and  1  detail levels. 
%  dwtr = [C; D]
%  C holds the coefficients (row vector)
%  D holds the details (row vectoe)
%
%   Example of Use:
%   data = [1; 0; -3; 2; 1; 0; 1; 2]; filter = [sqrt(2)/2 sqrt(2)/2];
%   wt = dwtr_1(data, filter)
%--------------------------------------------------------------------------

n = length(filterh);                  % Length of wavelet filter
C = data';                            % Data (row vector) live in V_j
dwtr = [];                            % At the beginning dwtr empty                
H  =  fliplr(filterh);                % Flip because of convolution
G  =  filterh;                        % Make quadrature mirror
G(1:2:n) = -G(1:2:n);                 %     counterpart
 
nn = length(C);                       % Length needed to make periodic 
C = [C(mod((-(n-1):-1),nn)+1)  C];    % Make periodic
D = conv(C,G);                        % Convolve,
D = D([n:2:(n+nn-2)]+1);              %     keep periodic and decimate
C = conv(C,H);                        % Convolve,
C = C([n:2:(n+nn-2)]+1);              %     keep periodic and decimate

dwtr = [C'; D'];
s = length(dwtr);
dwtr = [dwtr; data(s+1:end)'];

%======================================================================

function index = get_slice_index(n, i, k)
index = '';
for j=1:i-1,  
    if(j==1)
        index = ':';
    else
        index = [index ',:'];
    end
end
if(i==1)
    index = [index num2str(k)];
else
    index = [index ',' num2str(k)];
end
for j=i+1:n,
    index = [index ',:'];
end

%======================================================================
function index = get_space_index(ssize);

n = length(ssize);
index = '';
for j=1:n,  
    if(j~=1)
         index = [index, ','];
    end    
    index = [index '1:' num2str(ssize(j))];
end


