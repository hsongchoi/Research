close all;
clear all;
warning off;
tic
num_sims = 1000;
% NN = [512 1024 2048];
% SS = [3 5 7 10];
NN=1024;
SS=3
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
    wd = FWT_PO( data, coarsest, qmf );

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
%   tau=(max(var(data)-1/mu,10^(-6)))^0.5;
    tau=(max((var(data)-1/mu)/(2*(1-pi_0)^2),10^(-6)))^0.5;

    %-------------------------------
    
    sigS=zeros(1,N);
    sigS(sm_ind) = wd(sm_ind);
    %------------------------------
    Bs=b_factor(wd(a:N),mu,tau);
    [bsorted bindex] = sort(Bs);
    
    SigBest=sigS;
    sigS=SigBest;
    
    for j = 1:length(Bs)
        sigS(a+bindex(j)-1) = wd(a+bindex(j)-1);
        reC = IWT_PO(sigS,coarsest ,qmf);
        mSe(j)=sum((signal-reC).^2)/N; 
    end;
    [minMse jind] = min(mSe);
    alpha(k) = bsorted(jind);
    Ba = bindex(1:jind);
    SigBest(a+Ba-1) = wd(a+Ba-1);
    
    ResSig=IWT_PO(SigBest,coarsest ,qmf);
    smodata = smodata + ResSig;
    
    MSqE=[MSqE, 1/N*sum((signal-ResSig).^2)];
    clear jind Ba bA mSe SigBest sigS;
end;
smodata=smodata/num_sims;
biassq =1/N * sum((signal - smodata).^2);
variance=mean(MSqE)-biassq;
mean=mean(MSqE)
%eval(['save C:\MATLAB\work\new_fdr\BFDRresults\blocks',num2str(N),'snr',num2str(SNR),'.mat']); 

toc
