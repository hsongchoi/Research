close all;
clear all;
warning off;
tic
N=512;
SNR=3;
alpha = 0.05;
pi_0 = 0.9;

if N==512 
    num_lev=6; 
elseif N==1024 
    num_lev=7; 
else 
    num_lev=8; 
end;

coarsest = log2(N) - num_lev;  

num_sims = 1000;
pi_1 = 1- pi_0;


MSqE1=[];
biassq1=[];
smodata1=zeros(1,N);
for k = 1:num_sims
    sig1 = zeros(1,N);
  
    
    %------data generation
    %signal1 = MakeSignal('Blocks', N);
        t = (1:N) ./N;
     pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
        hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
        signal1 = zeros(size(t));
        for j=1:length(pos)
            signal1 = signal1 + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
        end

    sigma1 = std(signal1);
    signal1 = signal1.* SNR/sigma1;
  
    noise1 = randn(1,N);
    data1 = signal1+noise1;

    %--------------Wavelet decomposition
    %qmf1 = MakeONFilter('Haar',4);
    qmf = [1 1] ./ sqrt(2);
    %-----------wd1 = dwtr(data1, num_lev, qmf)---------
    n = length(qmf);                 
    C = data1(:)';                          
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
    %---------------------------
    sig1(sm_ind) = wd1(sm_ind);
    %--------mu&tau--------------------
    finest_lev1=wd1(fin_ind);
    q1_1=prctile(finest_lev1,25);
    q2_1=prctile(finest_lev1,75);
    pseudos1 = abs(q2_1-q1_1)/1.5;  
    mu1 = 1/pseudos1^2;        
    aa1=3*(var(data1)-1/mu1);
    b=(1-pi_0)^2;
    m1=(max((aa1/b),10^(-6)))^0.5;
    h1=sqrt(2*mu1);
    %----------Bayesian Factor for wavelet coefs (all but smoothe);
    if (-m1 > wd1(a:N))
    num= h1/2 .* exp(-h1 .* abs(wd1(a:N)));
    denom=(exp(h1.*wd1(a:N)).*sinh(h1.*m1))./(2*m1) ;
    Bs1=num./denom;
    elseif  (wd1(a:N)>= -m1) & (m1 >= wd1(a:N))
    num= h1/2 .* exp(-h1 .* abs(wd1(a:N)));
    denom=1/(2*m1)-exp(-h1.*m1).*cosh(h1.*wd1(a:N))./(2*m1);
    Bs1=num./denom;
    else 
    num= h1/2 .* exp(-h1 .* abs(wd1(a:N)));
    denom=(exp(-h1.*wd1(a:N)).*sinh(h1.*m1))./(2*m1);
    Bs1=num./denom;
    end

    %------------------------
    Ps1 = Bs1./(pi_1/pi_0 + Bs1);
 
   %----------------------------------
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
    
    R1=R1-1;
    s_ind1=Ind_sort1(1:R1);
    sig1(s_ind1+(a-1))=wd1(s_ind1+(a-1));
    %ResSig1 = idwtr(sig1, num_lev, qmf);
    nn = length(sig1);   n = length(qmf);  
    %if nargin==2, num_lev = round(log2(nn)); end;        
      H = fliplr(qmf);                                   
      G = qmf; G(2:2:n) = -G(2:2:n);
    LL = nn/(2^num_lev);                                   
    C =  sig1(1:LL);                                
    for j = 1:num_lev                                      
       w  = mod(0:n/2-1,LL)+1;                      
       D  = sig1(LL+1:2*LL);                       
       Cu(1:2:2*LL+n) = [C C(1,w)];                   
       Du(1:2:2*LL+n) = [D D(1,w)];                    
       C  = filter(H,[1],Cu) + filter(G,[1],Du); 
       C  = C([n:n+2*LL-1]-1);                       
       LL = 2*LL;                                      
    end;
    ResSig1 = C;            
    smodata1 = smodata1 + ResSig1;
    MSqE1=[MSqE1, 1/N*sum((signal1-ResSig1).^2)];

end;

smodata1=smodata1/num_sims;
biassq1 =1/N * sum((signal1 - smodata1).^2);
variance1=mean(MSqE1)-biassq1;
mean1=mean(MSqE1)
