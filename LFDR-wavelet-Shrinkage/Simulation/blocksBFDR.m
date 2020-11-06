% GAMMA 2.5
close all;
clear all;
warning off;
tic
%all_data=load('c:\matlab6p5\work\fdr\afm.dat');

addpath('C:\Users\haesong\Desktop\Haesong\Multiple\fdr')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Orthogonal')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Utilities')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Datasets')
addpath('C:\Users\haesong\Desktop\Haesong\Multiple\new_fd_matlabr')
num_sims = 1000;
%NN = [512 1024 2048];
NN=1024;
SS=7
%SS = [3 5 7 10];

     N=NN(1);
    SNR=SS(1);
if N==512 
    num_lev=4; 
elseif N==1024 
    num_lev=5; 
else 
    num_lev=6; 
end;
gam = 2.5;
coarsest = log2(N) - num_lev;
finest_level = log2(N) - 1;
MSqE=[];
biassq=[];
smodata=zeros(1,N);
for k=1:num_sims
    %------------data generation
    signal = MakeSignal('HeaviSine', N);
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
    %--------mu----------!!!!!!!!!!!!!!!!---------
    finest_lev=wd(fin_ind);

    q1=prctile(finest_lev,25);
    q2=prctile(finest_lev,75);
    pseudos = abs(q2-q1)/1.5; 
    mu = 1/pseudos^2;      
    %----------------
    
    sigS=zeros(1,N);
    sigS(sm_ind) = wd(sm_ind);
    %------------------------------

   j = 1; 
   for i = finest_level:-1:coarsest
    
    PI(i) = 1 - 1/(i-coarsest +1).^gam;
   tau(i)=(max((var(data)-1/mu)/(2*(1-PI(i))^2),10^(-6)))^0.5;
    %tau=(max(var(data)-1/mu,10^(-6)))^0.5;
    lev_ind=dyad(log2(N)-j);
    Bs = bb_factor(wd(lev_ind), mu, tau(i));
    bi = find(Bs < 1);
    sigS(lev_ind(1)+ bi -1) =wd(lev_ind(1)+bi - 1);  
    j = j+1;
   end;
    %-------------------------------
       
    ResSig=IWT_PO(sigS, coarsest,qmf);
    smodata = smodata + ResSig;
    
    MSqE=[MSqE, 1/N*sum((signal-ResSig).^2)];
end;
smodata=smodata/num_sims;
biassq =1/N * sum((signal - smodata).^2);
variance=mean(MSqE)-biassq;
mean=mean(MSqE)
%eval(['save C:\MATLAB\work\new_fdr\BFDRresults\blocks_2_',num2str(N),'snr',num2str(SNR),'.mat']); 

toc
