function g=dualshg(x)
% Usage:
% g=dualshg(x) 
% Returns the gradient <g> of the exact penalty function 
% at a point <x> for Shell Dual problem
global A B C D E Penalty_Coeff
if size(x,2)>1, x=x';end
g=x; pen=max([1.e3,Penalty_Coeff]);
g(1:5)=6*D.*x(1:5).^2 + 2*C*x(1:5);
g(6:15)= - B;
f=0;k=0;
for i=1:5
 d=A(:,i)'*x(6:15) - 2*C(i,:)*x(1:5) - 3*D(i)*x(i)^2 - E(i);
 if d>f, f=d; k=i; end
end
for i=1:15
 if -x(i)>f, f=-x(i); k=i+5; end
end
if k>0
 if k>5, g(k-5)=g(k-5)-pen;
 else, g(6:15)=g(6:15)+pen*A(:,k);
       g(1:5)=g(1:5)-pen*2*C(:,k);
       g(k)=g(k)-pen*6*D(k)*x(k);
 end
end
