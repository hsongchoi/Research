%test
close all
clear all
all_data=load('C:\Users\haesong\Desktop\Haesong\Multiple\afm.dat');
addpath('C:\Users\haesong\Desktop\Haesong\Multiple\fdr')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Orthogonal')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Utilities')
addpath('C:\Users\haesong\Desktop\WAVELAB850\Wavelab850\Datasets')
all_data=load('C:\Users\haesong\Desktop\Haesong\Multiple\afm.dat');
randn('seed',1)
lev = 10;
k = 2;
sizen = 1.;
%fi = MakeONFilter('Symmlet',4);
%si = MakeSignal('Doppler', 2^lev);
fi = MakeONFilter('Symmlet',4);
si = MakeSignal('HeaviSine', 2^lev);
sino = si + sizen * randn(size(si));
xx = linspace(0,1, length(si));
figure(1), plot(xx, sino)
%--------------------------------
siwa= FWT_PO(sino, k, fi);
smooth = siwa(1:2^k);
details =  siwa(2^k+1:2^lev);
%--------------------------------
finest = siwa(dyad(lev-1));
sighat = std(finest);
%--------------------------------
pvals=[];
nn = length(details);
dmin = mean(details);
for i = 1:nn
   pvals = [pvals 2*(1-cdf_nor(abs( (details(i)-dmin)/sighat ))) ];
end
figure(2)
hist(pvals, 30)
%--------------------------------
%tippet - Willcoxon Method
pvalss = sort(pvals);
LL=[];
for i = 1:nn
    LL = [LL  1-cdf_f( i/(nn-i+1) * (1- pvalss(i))/pvalss(i), 2*(nn-i+1), 2*i) ];
end
figure(3)
plot(LL)
    inde = min(find(LL >= 0.01))
    critp = pvalss(inde)
    mask = pvals < critp;
    thresh1 = [smooth mask .* details];
    sigt = IWT_PO(thresh1, k, fi);
figure(4)
plot(xx, sigt)
%---------------------------------
% Liptak Stoufer
%
LLL =[];
maska = [];
for i=1:nn
    shouldbenormal = perc_nor(pvals(i)+0.0001,0,1);
    LLL =[LLL shouldbenormal];
    maska = [maska abs(shouldbenormal) > sqrt(1.5*log(nn))];
end
 thresh2 = [smooth  maska .* details];
    sigt2 = IWT_PO(thresh2, k, fi);
figure(5)
plot(xx, sigt2)
hold on
plot(xx, si,'r-')