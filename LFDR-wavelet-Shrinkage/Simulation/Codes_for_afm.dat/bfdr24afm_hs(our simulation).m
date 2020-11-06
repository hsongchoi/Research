clear all;
close all;

%all_data=load('c:\matlab6p5\work\fdr\afm.dat');
%all_data=load('c:\matlab6p5\work\fdr\afm.dat');
all_data=load('C:\Users\haesong\Desktop\Haesong_Brani\fdr\afm.dat');
addpath('C:\Users\haesong\Desktop\Haesong_Brani\fdr')
addpath('C:\Users\haesong\Desktop\Haesong_Brani\WAVELAB850\Wavelab850\Orthogonal')
addpath('C:\Users\haesong\Desktop\Haesong_Brani\WAVELAB850\Wavelab850\Utilities')


sig=all_data(1:2048)'; clear all_data;
alpha=0.05;
N=2048;
num_lev=4

qmf = MakeONFilter('Symmlet', 4);

wd = FWT_PO( sig, num_lev, qmf );
varsig= var(sig)
%---------------------------

for j = 1:10
	eval(sprintf('fin_%d = dyad(log2(N)-%d);', j, j));
    eval(sprintf('wd_%d = wd(fin_%d);', j, j)); 
    eval(sprintf('std_%d = std(wd_%d);', j, j)); 
    eval(sprintf('var(%d) = std_%d^2;', j, j)); %OH
   eval(sprintf('clear std_%d  ;', j));
   eval(sprintf('dif =var-varsig;', j)); 
end
dif %8th=max 0.0174
SNR=dif(8)/varsig
m=sqrt(3*SNR)/(1-0.95) %251.3018


finest_lev=wd(fin_1); % generate the 1025th~2048th wd

q1=prctile(finest_lev,25);
q2=prctile(finest_lev,75);
pseudos = abs(q2-q1)/1.5; %sd
mu = 1/pseudos^2;        

qq1=prctile(wd, 25);
qq2=prctile(wd, 75);
ppseudos=abs(qq2-qq1)/1.5;
sigmatheta=(ppseudos^2)-(pseudos^2);
%tau=sqrt((N-2)/N)*sigmatheta;
Bs=b_factor(wd,mu,m); 

Ps=Bs./(0.01+0.99.*Bs);  


[P_sort, Ind_sort]=sort(Ps); %P_sort: the sorted Ps, Ind_sort: Index of the sorted Ps
EQ=0;R=1;
while EQ<alpha
    R=R+1;
    P0ly=[1-P_sort(1), P_sort(1)];
    for i=2:R
        b=[1-P_sort(i), P_sort(i)];
        P0ly=conv(P0ly,b);  
    end;
    P0ly1=P0ly(2:R+1);
    zz=1:R;
    P0ly1=zz.*P0ly1;
    EQ=sum(P0ly1)/R;
    EQQ(R)=EQ;
end;
R=R-1;
s_ind=Ind_sort(1:R);
smsig=zeros(1,N);
smsig(s_ind)=wd(s_ind);
ResSig=IWT_PO(smsig,num_lev,qmf);
bs=[]


plot(sig,'c'); hold on;  
plot(ResSig)