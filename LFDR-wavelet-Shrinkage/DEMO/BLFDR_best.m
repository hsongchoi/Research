close all;
clear all;
warning off;
tic
num_sims = 1000;
N=512;
SNR=3;
if N==512 
    num_lev=4; 
elseif N==1024 
    num_lev=5; 
else 
    num_lev=6; 
end;
coarsest = log2(N) - num_lev;
MSqE=[];
biassq=[];
smodata=zeros(1,N);
pi_0 = 0.95;
for k=1:num_sims
    %------------data generation
    %signal = MakeSignal('Blocks', N);
    t = (1:N) ./N;
     pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
        hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
        signal = zeros(size(t));
        for j=1:length(pos)
            signal = signal + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
        end

    sigma = std(signal);
    signal = signal.*SNR/sigma;
    noise = randn(1,N);
    data = signal + noise;
    
    %-------------Wavelet Decomposition
     
    %qmf = MakeONFilter('Haar', 4);
    qmf = [1 1] ./ sqrt(2);
    %-----------wd1 = dwtr(data1, num_lev, qmf)---------
    n = length(qmf);                 
    C = data(:)';                          
    dwtr = [];                                                 
    H = qmf;
    G = fliplr(qmf);
    G(1:2:n) = -G(1:2:n);             
    for j = 1:num_lev                             
       nn = length(C);                    
       C = [C(mod((-(n-1):-1),nn)+1)  C]; 
       D=filter(G,[1],C);                 

       D = D([n:2:(n+nn-2)]+1);   
       C=filter(H,[1],C);               
       C = C([n:2:(n+nn-2)]+1);     
       dwtr = [D,dwtr];                                        
    end;                                      
       wd = [C, dwtr];   
    %---------------------------
    fin_ind = (2^(log2(N)-1)+1):(2^(log2(N))) ;
    a = (2^(log2(N)-num_lev)+1):(2^(log2(N)-num_lev+1)) ;
    sm_ind=1:a(1)-1; a=a(1);
    %--------mu&tau-----------!!!!!!!!!!!!!!!!---------
    finest_lev=wd(fin_ind);

    q1=prctile(finest_lev,25);
    q2=prctile(finest_lev,75);
    pseudos = abs(q2-q1)/1.5; 
    mu = 1/pseudos^2;      
    aa=3*(var(data)-1/mu);
    b=(1-pi_0)^2;
    m=(max((aa/b),10^(-6)))^0.5;
    h=sqrt(2*mu);

    %-------------------------------
    
    sigS=zeros(1,N);
    sigS(sm_ind) = wd(sm_ind);
    %------------------------------
    %Bs=b_factor(wd(a:N),h,m);
    if (-m > wd(a:N))
    num= h/2 .* exp(-h.* abs(wd(a:N)));
    denom=(exp(h.*wd(a:N)).*sinh(h.*m))./(2*m) ;
    Bs=num./denom;
    elseif  (wd(a:N)>= -m) & (m >= wd(a:N))
    num= h/2 .* exp(-h .* abs(wd(a:N)));
    denom=1/(2*m)-exp(-h.*m).*cosh(h.*wd(a:N))./(2*m);
    Bs=num./denom;
    else 
    num= h/2 .* exp(-h .* abs(wd(a:N)));
    denom=(exp(-h.*wd(a:N)).*sinh(h.*m))./(2*m);
    Bs=num./denom;
    end

    [bsorted bindex] = sort(Bs);
    
    SigBest=sigS;
    sigS=SigBest;
    
    for j = 1:length(Bs)
        sigS(a+bindex(j)-1) = wd(a+bindex(j)-1);
        %reC = idwtr(sigS, num_lev, qmf);
        nn = length(sigS);   n = length(qmf);  
        %if nargin==2, num_lev = round(log2(nn)); end;        
          H = fliplr(qmf);                                   
          G = qmf; G(2:2:n) = -G(2:2:n);
        LL = nn/(2^num_lev);                                   
        C =  sigS(1:LL);                                
        for j = 1:num_lev                                      
       w  = mod(0:n/2-1,LL)+1;                      
       D  = sigS(LL+1:2*LL);                       
       Cu(1:2:2*LL+n) = [C C(1,w)];                   
       Du(1:2:2*LL+n) = [D D(1,w)];                    
       C  = filter(H,[1],Cu) + filter(G,[1],Du); 
       C  = C([n:n+2*LL-1]-1);                       
       LL = 2*LL;                                      
    end;
    reC = C;       
        mSe(j)=sum((signal-reC).^2)/N; 
    end;
    [minMse jind] = min(mSe);
    alpha(k) = bsorted(jind);
    Ba = bindex(1:jind);
    SigBest(a+Ba-1) = wd(a+Ba-1);
    
    %ResSig = idwtr(SigBest, num_lev, qmf);
    nn = length(SigBest);   n = length(qmf);
    %if nargin==2, num_lev = round(log2(nn)) end;        
      H = fliplr(qmf);                                   
      G = qmf; G(2:2:n) = -G(2:2:n);
    LL = nn/(2^num_lev);                                   
    C =  SigBest(1:LL);                                
    for j = 1:num_lev                                      
       w  = mod(0:n/2-1,LL)+1;                      
       D  = SigBest(LL+1:2*LL);                       
       Cu(1:2:2*LL+n) = [C C(1,w)];                   
       Du(1:2:2*LL+n) = [D D(1,w)];                    
       C  = filter(H,[1],Cu) + filter(G,[1],Du); 
       C  = C([n:n+2*LL-1]-1);                       
       LL = 2*LL;                                      
    end;
    ResSig = C;       
    smodata = smodata + ResSig;
    
    MSqE=[MSqE, 1/N*sum((signal-ResSig).^2)];
    clear jind Ba bA mSe sigS;
end;
smodata=smodata/num_sims;
biassq =1/N * sum((signal - smodata).^2);
variance=mean(MSqE)-biassq;
mean=mean(MSqE)
