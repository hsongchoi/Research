close all;
clear all;
warning off;
tic
num_sims = 1000;
NN=512;
SS=3
%NN = [512 1024 2048];
%SS = [3 5 7 10];

%addpath('C:\Users\haesong\Desktop\Haesong\Multiple\fdr')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Orthogonal')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Utilities')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Datasets')
addpath('C:\Users\haesong\Desktop\Haesong\Multiple\new_fd_matlabr')
%for t=1:3
 %   N=NN(t);
%for tt=1:4
 %   SNR=SS(tt);
    N=NN(1);
    SNR=SS(1);
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
    signal = MakeSignal('Blocks', N);
    sigma = std(signal);
    signal = signal.*SNR/sigma;
    noise = randn(1,N);
    data = signal + noise;
    
    %-------------Wavelet Decomposition
     
    qmf = MakeONFilter('Haar', 4);
    wd = FWT_PO( data, coarsest , qmf );

    %---------------------------
    fin_ind=dyad(log2(N)-1);
    a=dyad(log2(N)-num_lev);
    sm_ind=1:a(1)-1; a=a(1);
    %--------mu&tau-----------!!!!!!!!!!!!!!!!---------
    finest_lev=wd(fin_ind);

    q1=prctile(finest_lev,25);
    q2=prctile(finest_lev,75);
    pseudos = abs(q2-q1)/1.5; 
    mu = 1/pseudos^2;      
   % m=sqrt(3*SNR)/(1-pi_0)
    m= 2*sqrt(SNR)
     %m=50;
    %-------------------------------
    
    sigS=zeros(1,N);
    sigS(sm_ind) = wd(sm_ind);
    %------------------------------
    wd(a:N)
    %Bs=b_factor(wd(a:N),mu,m);
    d=wd(a:N);
    h=mu
    if (d <= -m)
    num= h/2 .* exp(-h .* abs(d));
denom=(exp(h.*d).*sinh(h.*m))./(2*m) ;
 output=num./denom;
elseif (-m <= d) & (d <= m)
       num= h/2 .* exp(-h .* abs(d));
denom=1/(2*m)-exp(-h.*m).*cosh(h.*d)./(2*m);
 output=num./denom;
elseif (m <= d)
     num= h/2 .* exp(-h .* abs(d));
denom=(exp(-h.*d).*sinh(h.*m))./(2*m);
 output=num./denom;
end
    
    Ba = find(Bs < 1);
    sigS(a+Ba-1) = wd(a+Ba-1);
    
    ResSig=IWT_PO(sigS, coarsest,qmf);
    smodata = smodata + ResSig;
    
    MSqE=[MSqE, 1/N*sum((signal-ResSig).^2)];
    clear jind Ba bA mSe SigBest sigS;
end;
smodata=smodata/num_sims;
biassq =1/N * sum((signal - smodata).^2);
variance=mean(MSqE)-biassq;
mean=mean(MSqE)
%eval(['save C:\MATLAB\work\new_fdr\BFDRresults\blocks_fixed',num2str(N),'snr',num2str(SNR),'.mat']); 

toc
