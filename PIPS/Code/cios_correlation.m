%%correlation
%c1:Course: How prepared to take subject
%c2:Course: Amount learned
%c3:Course: Assignments facilitated learning
%c4:Course: Assignments measured knowledge
%c5:Course: Overall effectiveness
%i1:Instructor: Clarity /c6
%i2:Instructor: Communicated how to succeed/c7
%i3:Instructor: Respect for students/c8
%i4:Instructor: Enthusiasm/c9
%i5:Instructor: Stimulates interest/c10
%i6:Instructor: Availability/c11
%i7:Instructor: Feedback helpfulness/c12
%i8:Instructor: Overall effectiveness/c13
%s1:Structure to work together
%s2:Structured for feedback
%s3:Structured for active students
clear all;
close all;
addpath('C:\Users\haesong\Desktop\8_1\cios');
data=load('C:\Users\haesong\Desktop\8_1\cios\comp3.csv');
c1=data(:,1);
c2=data(:,2);
c3=data(:,3);
c4=data(:,4);
c5=data(:,5);
c6=data(:,6);
c7=data(:,7);
c8=data(:,8);
c9=data(:,9);
c10=data(:,10);
c11=data(:,11);
c12=data(:,12);
c13=data(:,13);
s1=data(:,14);
s2=data(:,15);
s3=data(:,16);

for j = 1:13
    	eval(sprintf('[r%d, p%d] = corr(s1,c%d);', j, j, j));
   
end
        r11=[r1, r2, r3, r4,r5, r6, r7,r8,r9,r10,r11,r12,r13]
        p11=[p1, p2, p3, p4, p5,p6,p7,p8,p9,p10,p11,p12,p13]
        
for j = 1:13
    	eval(sprintf('[r%d, p%d] = corr(s2,c%d);', j, j, j));
   
end
 r22=[r1, r2, r3, r4,r5, r6, r7,r8,r9,r10,r11,r12,r13]
        p22=[p1, p2, p3, p4, p5,p6,p7,p8,p9,p10,p11,p12,p13]
        
 for j = 1:13
    	eval(sprintf('[r%d, p%d] = corr(s3,c%d);', j, j, j));
   
 end
 r33=[r1, r2, r3, r4,r5, r6, r7,r8,r9,r10,r11,r12,r13]
        p33=[p1, p2, p3, p4, p5,p6,p7,p8,p9,p10,p11,p12,p13]
        
 tbl1=table(s1,s2,s3, c5) %c5:Course: Overall effectiveness
 lm1 = fitlm(tbl1,'c5~s1+s2+s3')
 plotResiduals(lm1,'probability')
 %%
 the linear regression analysis requires all variables to be multivariate normal.
 [s11,lambda] = boxcox(s1);
  [s22,lambda] = boxcox(s2);
    [s33,lambda] = boxcox(s3);
      [c55,lambda1] = boxcox(c5);
      [c1313, lambda2]=boxcox(c13);
      
      imp=(c5.^lambda1-1)./lambda1 %cox-transformation
      
      
    tbl11=table(s1,s2,s3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c55,c1313);
       lm2 = fitlm(tbl11,'c55~s1+s2+s3+c1+c2+c3+c4+c6+c7+c8+c9+c11+c12')
      lm1 = fitlm(tbl11,'c5~s1+s2+s3+c1+c2+c3+c4+c6+c7+c8+c9+c11+c12')
      lm3= fitlm(tbl11,'c13~s1+s2+s3+c1+c2+c3+c4+c6+c7+c8+c9+c11+c12')
      lm4= fitlm(tbl11,'c1313~s1+s2+s3+c1+c2+c3+c4+c6+c7+c8+c9+c11+c12')
       %i1:Instructor: Clarity /c6
%i2:Instructor: Communicated how to succeed/c7
%i3:Instructor: Respect for students/c8
%i4:Instructor: Enthusiasm/c9
%i5:Instructor: Stimulates interest/c10
%i6:Instructor: Availability/c11
%i7:Instructor: Feedback helpfulness/c12
%i8:Instructor: Overall effectiveness/c13
%s1:Structure to work together
%s2:Structured for feedback
%s3:Structured for active students
plotResiduals(lm1,'probability')
       plotResiduals(lm2,'probability')
        plotResiduals(lm3,'probability') 
          plotResiduals(lm4,'probability')
       plotResiduals(lm2)
       kstest(lm1.Residuals)
       The probability plot that results is more linear
       xy=[s1, s2, s3 c55]
       VIF = SingleVIF(lm2)
       x=[s1, s2 ,s3]
       [ndim, prob] = barttest(x,0.05) %Bartlett’s test. Tests if the variances of the data values along each principal component are equal, against the alternative that the variances are not all equal.


addpath('D:\cios')
lm= regstats(c55, [s1 s2 s3], 'linear', {'r','yhat'}) 
TestHet(lm.r,[s1, s2 s3], '-Ws', lm.yhat)
%%
sqrt=sqrt(c5)
lm= regstats(sqrt, [s1 s2 s3], 'linear', {'r','yhat'}) 
TestHet(lm.r,[s1, s2 s3], '-Ws', lm.yhat)

log=log(c5)
lm2= regstats(log, [s1 s2 s3], 'linear', {'r','yhat'}) 
TestHet(lm2.r,[s1, s2 s3], '-BPK')

re=asin(c5.^(1/2))
lm2= regstats(re, [s1 s2 s3], 'linear', {'r','yhat'}) 
TestHet(lm2.r,[s1, s2 s3], '-BPK')
[b,dev,stats] = glmfit([s1 s2 s3],c55,'link','log','weights')
