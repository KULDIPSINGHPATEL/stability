
%% symbolic computation of BetaN'
function [D1N]=betaandD_new(N)
syms(sym('A',[1 N-1]));A=(sym('A',[1 N-1]));syms(sym('B',[1 N-1]));
B=(sym('B',[1 N-1]));syms(sym('C',[1 N-1]));C=(sym('C',[1 N-1]));
syms(sym('P',[1 N+1]));P=(sym('P',[1 N+1]));P(1)=1;
  P(3)=-B(1)*P(2)-C(1);
  j=4;
  for i=1:N-2
  P(j)=(-C(j-2))*A(j-3)*P(j-2)+(-B(j-2))*P(j-1);
  j=j+1;
  end
 %BetaN=P(N+1); % the expression of product of A_js and BetaN'
 K=coeffs(P(N+1),P2);
 D1N=simplifyFraction(K(2)); % expresssion for D1N
end