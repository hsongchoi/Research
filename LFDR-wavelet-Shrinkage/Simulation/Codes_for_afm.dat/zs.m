%z****************************'s
clear all;
close all;
%--------------
sizn = 0.1;
deep = 3;
nn = 1024;
%-----------------------------
sig = MakeSignal('Doppler', nn);
sino = sig + sizn * randn( size(sig));
%----------------------------------
filt = MakeONFilter('Symmlet', 4);
wd = FWT_PO( sino, deep, filt );
   stdd = std( wd(dyad(log2(nn)-1)));
%----------------------------------
pvals = 2 .* (1 - cdf_nor( abs( (wd-mean(wd))/stdd )));
zs = perc_nor(pvals+0.000001, 0, 1);
figure(1)
hist(zs, 30)
%---------------------------------
[aa, bb] = hist(zs, 100);
 h1=0.15;
 h2=0.30;
 h3=0.60;
 xx=sort(bb);
 figure(2)
 plot(xx,aa,'*')
 lle1=[];
 lle2=[];
 lle3=[];
 for j=-480:310
  lle1=[lle1 loc_lin(j/100, xx, aa, h1)];
  lle2=[lle2 loc_lin(j/100, xx, aa, h2)];
  lle3=[lle3 loc_lin(j/100, xx, aa, h3)];
 end
 hold on
 plot((-480:310)/100, lle1, 'r--' )
 plot((-480:310)/100, lle2, 'k-.' )
 plot((-480:310)/100, lle3, 'b-')
 legend('data','h=0.15','h=0.30','h=0.60',1)
 hold off
 %---------------------------------------------
 xxx = xx( xx > -3.3);
 aaa = aa( xx > -3.3);
 aaa = length(xx)/length(xxx) * aaa;
  figure(3)
 plot(xxx,aaa,'*')
 lle1=[];
 lle2=[];
 lle3=[];
 for j=-330:310
  lle1=[lle1 loc_lin(j/100, xxx, aaa, h1)];
  lle2=[lle2 loc_lin(j/100, xxx, aaa, h2)];
  lle3=[lle3 loc_lin(j/100, xxx, aaa, h3)];
 end
 hold on
 plot((-330:310)/100, lle1, 'r--' )
 plot((-330:310)/100, lle2, 'k-.' )
 plot((-330:310)/100, lle3, 'b-')
 %-------------
 ss = (-330:310)/100;
 dd = 1/sqrt(2 * pi) * exp( - 1/2 * ss.^2);
 plot((-330:310)/100, 128 * dd, 'b-','linewidth', 3)
 legend('data','h=0.15','h=0.30','h=0.60',1)
 hold off
 %---------------------------------------------
