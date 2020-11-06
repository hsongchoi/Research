% GAMMA 2
close all;
clear all;
warning off;
tic

num_sims = 1000;
alpha=0.05;
N=1024;
SNR=7;
gam = 2;
phi_0 = 0.95;
phi_1 = 1- phi_0;
if N==512 
    num_lev=4; 
elseif N==1024 
    num_lev=5; 
else 
    num_lev=6; 
end;

coarsest = log2(N) - num_lev;
finest_level = log2(N) - 1;
MSqE=[];
biassq=[];
smodata=zeros(1,N);
for k=1:num_sims
    %------------data generation
    %signal1 = MakeSignal('Blocks', N);
            t = (1:N) ./N;
     pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
        hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
        signal1 = zeros(size(t));
        for j=1:length(pos)
            signal1 = signal1 + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
        end

    sigma = std(signal1);
    signal1 = signal1.*SNR/sigma;
    noise = randn(1,N);
    data = signal1 + noise;
    
    %-------------Wavelet Decomposition
     
    %qmf = MakeONFilter('Haar', 4);
    qmf = [1 1] ./ sqrt(2);
    %-----------wd1 = dwtr(data, num_lev, qmf)---------
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
       wd1 = [C, dwtr];   
    %---------------------------
    fin_ind = (2^(log2(N)-1)+1):(2^(log2(N))) ;
    a = (2^(log2(N)-num_lev)+1):(2^(log2(N)-num_lev+1)) ;
    sm_ind=1:a(1)-1; a=a(1);
    %--------mu----------!!!!!!!!!!!!!!!!---------
    finest_lev=wd1(fin_ind);
    q1=prctile(finest_lev,25);
    q2=prctile(finest_lev,75);
    pseudos = abs(q2-q1)/1.5; 
    mu = 1/pseudos^2;      
    %----------------
    
    sigS=zeros(1,N);
    sigS(sm_ind) = wd1(sm_ind);
    %------------------------------

    j = 1; 
    for i = finest_level:-1:coarsest
    PI(i) = 1 - 1/(i-coarsest +1).^gam;
    aa=3*(var(data)-1/mu);
    b(i)=(1-PI(i))^2;
    m1(i)=(max((aa/b(i)),10^(-6)))^0.5;
    h1=sqrt(2*mu);
    lev_ind=dyad(log2(N)-j);
    %Bs=b_factor(wd1(lev_ind),h1,m1(i));
    if (-m1(i) > wd1(lev_ind))
    num= h1/2 .* exp(-h1 .* abs(wd1(lev_ind)));
    denom=(exp(h1.*wd1(lev_ind)).*sinh(h1.*m1(i)))./(2*m1(i)) ;
    Bs=num./denom;
    elseif  (wd1(lev_ind)>= -m1(i)) & (m1(i) >= wd1(lev_ind))
    num= h1/2 .* exp(-h1 .* abs(wd1(lev_ind)));
    denom=1/(2*m1(i))-exp(-h1.*m1(i)).*cosh(h1.*wd1(lev_ind))./(2*m1(i));
    Bs=num./denom;
    else 
    num= h1/2 .* exp(-h1 .* abs(wd1(lev_ind)));
    denom=(exp(-h1.*wd1(lev_ind)).*sinh(h1.*m1(i)))./(2*m1(i));
    Bs=num./denom;
    end

    Ps1 = Bs./(phi_1/phi_0 + Bs);
    %%%%%%%
    [P_sort1, Ind_sort1]=sort(Ps1);
    EQ1=0;R1=1;

    while EQ1<alpha
        R1=R1+1;
        P0ly_1=[1-P_sort1(1), P_sort1(1)];
        for i=2:R1
            b=[1-P_sort1(i), P_sort1(i)];
            P0ly_1=conv(P0ly_1,b);
        end;
        P0ly1_1=P0ly_1(2:R1+1);
        zz1=1:R1;
        P0ly1_1=zz1.*P0ly1_1;
        EQ1=sum(P0ly1_1)/R1;
        EQQ(R1)=EQ1;
    end;
   pmax=P_sort1(R1);
   thres=((pmax/(1-pmax))*(phi_1/phi_0));
   %thres=((pmax/(1-pmax))*0.05);
    bi = find(Bs < 1);
    sigS(lev_ind(1)+ bi -1) =wd1(lev_ind(1)+bi - 1);  
    j = j+1;
   end;
    %-------------------------------
       
    %ResSig = idwtr(sigS, num_lev, qmf);
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
    ResSig = C;   
    smodata = smodata + ResSig;
    
    MSqE=[MSqE, 1/N*sum((signal1-ResSig).^2)];
end;
smodata=smodata/num_sims;
biassq =1/N * sum((signal1 - smodata).^2);
mean=mean(MSqE)
variance=mean-biassq;

toc
