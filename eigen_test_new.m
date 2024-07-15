
%% To test that roots of D1N are positive
clc;clear;N=input('Enter value of N=');
z_l=input('Enter value of left space boundary point z_l=');
z_r=input('Enter value of right space boundary point z_r=');
dv=input('Enter value of time step dv=');
syms a(z) b(z) w1(z) y1(z) w2(z) y2(z)
a(z)=input('Enter the convection coefficient a(z)=');
b(z)=input('Enter the convection coefficient b(z)=');
x=linspace(z_l,z_r,N+1);dz=x(2)-x(1);%%% discretizing space domain
%%% defining derivatives of a(z) and b(z)
w=a(x);y=b(x);w1(z)=diff(a(z));y1(z)=diff(b(z));w2(z)=diff(w1(z));y2(z)=diff(y1(z));
wd=w1(x);yd=y1(x);wdd=w2(x);ydd=y2(x);
%%%% defining p q r l m n as per manuscript
 gamma=zeros(N-1,1);xeta=zeros(N-1,1);alpha=zeros(N-1,1);
 p=zeros(N-1,1);q=zeros(N-1,1);r=zeros(N-1,1);
 l=zeros(N-1,1);m=zeros(N-1,1);n=zeros(N-1,1);
 for i=1:N-1 
 gamma(i)=(w(i+1)/y(i+1))+(2/y(i+1))*yd(i+1);
 xeta(i)=w(i+1)-((dz^2)/12)*((w(i+1)/y(i+1))*wd(i+1)+(2/y(i+1))*wd(i+1)*yd(i+1)-wdd(i+1));
 alpha(i)=y(i+1)+((dz^2)/12)*(((w(i+1)/y(i+1))*yd(i+1))-2*wd(i+1)+ydd(i+1)...
     -((2/y(i+1))*(yd(i+1)^2))+((w(i+1)^2)/y(i+1)));
 p(i)=((2+dz)*gamma(i))/(24*dv);
 q(i)=5/(6*dv);
 r(i)=((2-dz)*gamma(i))/(24*dv);
 l(i)=(-xeta(i)/(2*dz))-(alpha(i)/(dz^2));
 m(i)=(2*alpha(i))/(dz^2);
 n(i)=(xeta(i)/(2*dz))-(alpha(i)/(dz^2));
 end
 [D1N]=betaandD_new(N);
 %fprintf('Take Beta_1_Prime=P2, then expression for Beta_N_Prime is= %s. \n',BetaN);
 fprintf('The expression for D1_N is= %s .\n',D1N);
 %%%% calling function to check the values
[u,Pol]=roots_value_new(N,p,q,r,l,m,n);
disp("The polynomial in lambda is=")
Pol
disp("The roots of the D1_N=")
u
if min(real(u)) > 0
disp('Result: Since the roots are positive, the compact Scheme is stable for above parameter')
else
 disp('Result: Since the roots are NOT positive, Stability condition is not met for above parameter')
end