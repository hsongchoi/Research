%z****************************'s
clear all;
close all;


 disp('DJ TEST FUNCTIONS')
  lw = 2; 
  set(0, 'DefaultAxesFontSize', 15);
  fs = 15;
  msize = 6;
 


%--------------
%randn('seed',1);
sizn = 0.1;
deep = 3;
nn = 1024;
th=0.3; %FDR threshold value
sigma=0.1; %sigma value for the estimation of f(z)
%-----------------------------
sig = MakeSignal('Doppler', nn);
sino = sig + sizn * randn( size(sig));
xxs = linspace(0,1, nn);
%----------------------------------
filt = MakeONFilter('Symmlet', 4);
wd = FWT_PO( sino, deep, filt );
stdd = std( wd(dyad(log2(nn)-1)));
thu = sqrt( 2 * log(nn) ) * stdd
%----------------------------------
pvals = 2 .* (1 - cdf_nor( abs( (wd-mean(wd))/stdd )));
zs = perc_nor(pvals+0.000001, 0, 1);
figure(1)
histfit(zs, 30)
%---------------------------------

[Z, IS]=sort(zs);
I=find((-3.5<Z)&(Z<3.5));
Zsel=Z(I);
mu_e=mean(Zsel);
sigma_e=std(Zsel);

f0_th=normpdf(zs,0,1);
f0_emp=normpdf(zs,mu_e,sigma_e);
f=ker_sm(zs,sigma);
fdr_th=f0_th./f;
fdr_emp=f0_emp./f;


figure(2);
plot(Z,f(IS),'linewidth', lw ); 
hold on
plot(Z,f0_th(IS),'r:','linewidth', lw );
hold on
plot(Z,f0_emp(IS),'g-.', 'linewidth', lw);
title(['sigma=',num2str(sigma),'  for estimated f(z)']);
legend('Estimated f(z)','Theoretical f_0(z)','Empirical f_0(z)',2);


figure(3)
plot(Z,fdr_th(IS), 'linewidth', lw);
hold on
plot(Z,fdr_emp(IS),'r', 'linewidth', lw);
hold on
plot(Z,th,'.', 'markersize', msize)
legend('Theoretical','Empirical',['fdr=',num2str(th)], 2);
title('Local false discovery rate');

%------------------------------------restoration
ind_th=find(fdr_th<th);
ind_emp=find(fdr_emp<th);

sigwdth_th=zeros(1,length(wd));
sigwdth_emp=zeros(1,length(wd));

sigwdth_th(ind_th)=wd(ind_th);
sigwdth_emp(ind_emp)=wd(ind_emp);

sigRes_th=IWT_PO(sigwdth_th, deep,filt);
sigRes_emp=IWT_PO(sigwdth_emp, deep,filt);

figure(4)
subplot(2,2,1);
plot(xxs,sig, 'linewidth', lw); title('original signal')
axis tight
subplot(2,2,2)
plot(xxs, sigRes_th, 'linewidth', lw); title('Restored signal (th)')
axis tight
subplot(2,2,3);
plot(xxs, sino, 'linewidth', lw); title('Noisy signal')
axis tight
subplot(2,2,4)
plot(xxs, sigRes_emp, 'linewidth', lw); title('Restored signal (emp)')
axis tight

errt= sum((sig - sigRes_th).^2)
erre= sum((sig - sigRes_emp).^2)

figure(5)
wdthu = wd .* (abs(wd) > thu);
recu = IWT_PO( wdthu, deep, filt );
plot( xxs, recu, 'linewidth', lw); title('Restored signal (uni)')
axis tight
erru= sum((sig - recu).^2)


figure(6)
plot(xxs, sigRes_th, 'linewidth', lw); title('Restored signal (th)')
axis tight

figure(7)
plot(xxs, sigRes_emp, 'linewidth', lw); title('Restored signal (emp)')
axis tight