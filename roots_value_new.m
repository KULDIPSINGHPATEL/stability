%% function to compute roots of D1N symbolically
function[u, D1N]=roots_value_new(N,p,q,r,l,m,n)
syms lambda
syms(sym('A',[1 N-1]));A=(sym('A',[1 N-1]));
syms(sym('B',[1 N-1]));B=(sym('B',[1 N-1]));
syms(sym('C',[1 N-1]));C=(sym('C',[1 N-1]));
%%% Now find in terms of lambda
   for h=1:N-1
       A(h)=n(h)-lambda*r(h);B(h)=m(h)-lambda*q(h);C(h)=l(h)-lambda*p(h);
   end
  syms(sym('P',[1 N+1]));P=(sym('P',[1 N+1]));P(1)=1;
  P(3)=-B(1)*P(2)-C(1);
  j=4;
  for i=1:N-2
  P(j)=(-C(j-2))*A(j-3)*P(j-2)+(-B(j-2))*P(j-1);
  j=j+1;
  end
 K=coeffs(P(N+1),P2); % K(1)= cosntant term and K(2)= Beta_1_prime
 D1N=simplifyFraction(K(2)); % expresssion for D1N (K(2) is taken to avoid the constant term)
 %%%% find coefficients of lambda
 M=coeffs(D1N,lambda); % first entry of M is zeroth order term and 2nd is coefifcient of lambda and so on
 %%% Arrange in reverse order to apply built in MATLAB function roots
 V=zeros(N,1);
 for b=1:N
     V(b)=M(N-b+1);
 end
 u=roots(V); %%%%%% find roots of D1N
end