function data = idwtrn(dwtr, L, filterh)
%  function data = idwtr_n(dwtr, L, filterh); 
%  Calculates the IDWT of wavelet transformation dwtr (n dimensional)
%  using scaling filter and L levels.
%--------------------------------------------------------------------------

data = dwtr;
d = ndims(dwtr);                                   % number of dimensions
bsize = size(dwtr);                                % size of the space

oned = ( d==2 && ~isempty(find(bsize==1)) );       % true, if vector (row or column)
 
if(oned)  msize = max(bsize);
else      msize= min(bsize);   end;

L = min([L floor(log2(msize))]);
for i = L:-1:1,                                    % for each level       
    ssize = floor(bsize/(2^(i-1)));                % new subspace size
    if(oned) ssize(find(ssize==0)) = 1; end
    sindex = get_space_index(ssize);               % subspace index 
    eval(['data(' sindex ') = idwtr_n1(data(' sindex '), filterh);']);    % replace subspace
end    


function data = idwtr_n1(dwtr, filterh)
% function data = idwtr_n1(dwtr, filterh); 
% Calculates the IDWT of wavelet transformation dwtr (n dimensional)
% using scaling filter  "filterh"  and  1 levels.  
%--------------------------------------------------------------------------

d = ndims(dwtr);                                   % dimension of the data
s = size(dwtr);                                    % size of each dimension

oned = ( d==2 && ~isempty(find(s==1)) );           % true, if vector (row or column)

if ( oned )                                        % one dimensional data, end of recursion            
    C_r = reshape(dwtr, [length(dwtr) 1]);         % reshape to one dimensional column vector
    C_t = idwtr_11(C_r, filterh);                  % call one one dimensional specific function    
    data = reshape(C_t, size(dwtr));               % reshape back to original shape
    return;
end

data = dwtr;                             
for i=d:-1:1,                                                 % for each dimension
    for j=1:s(i),                                             % for each slice within dimension 
        C_s = eval(['data(' get_slice_index(d, i, j) ')']);   % get the slice
        if (d==2)                                             % if the slice is one dimensional 
            C_r = reshape(C_s, [length(C_s) 1]);              % reshape to column vector    
        else                                                   
            C_r = squeeze(C_s);                               % remove the virtual dimension
        end
        C_t = idwtr_n1(C_r , filterh);                        % solve the n-1 dimensional problem 
        C_r = reshape(C_t, size(C_s));                        % reshape back to slice dimensions
        eval(['data(' get_slice_index(d, i, j) ') = C_r;']);  % set old slice to new slice  
    end
end



function  data = idwtr_11(dwtr, filterh)
% function data = idwtr_11(dwtr, filterh); Calculates the IDWT of wavelet
% transformation dwtr using scaling filter  "filterh"  and  1 levels.  
% dwtr is a column vector
%--------------------------------------------------------------------------
nn = length(dwtr);   n = length(filterh);       % Lengths
LL = floor(nn/2);                               % Number of scaling coeffs
H = filterh;                                    % Wavelet H filter
G = fliplr(H); G(2:2:n) = -G(2:2:n);            % Wavelet G filter

C =  dwtr(1:LL)';                               % Scaling coeffs

w  = mod(0:n/2-1,LL)+1;                         % Make periodic
D  = dwtr(LL+1:2*LL)';                          % Wavelet coeffs
Cu(1:2:2*LL+n) = [C C(1,w)];                    % Upsample & keep periodic
Du(1:2:2*LL+n) = [D D(1,w)];                    % Upsample & keep periodic
C  = conv(Cu,H) + conv(Du,G);                   % Convolve & add
C  = C([n:n+2*LL-1]-1);                         % Periodic part

data = [C'; dwtr(2*LL+1:end)];                     % The inverse DWT


%========================================================================
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

%========================================================================
function index = get_space_index(ssize);

n = length(ssize);
index = '';
for j=1:n,  
    if(j~=1)
         index = [index, ','];
    end    
    index = [index '1:' num2str(ssize(j))];
end



