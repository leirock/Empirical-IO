function f=dualshf(x)
% Usage:
% f=dualshf(x) 
% Returns the value of the exact penalty function <f>
% at a point <x> for the Shell Dual problem.
global A B C D E Penalty_Coeff
pen=max([1.e3,Penalty_Coeff]);
if size(x,2)>1, x=x';end
f0=2*D'*x(1:5).^3 + (C*x(1:5))'*x(1:5) - B'*x(6:15);
f=max([0;A'*x(6:15) - 2*C*x(1:5) - 3*D.*x(1:5).^2 - E; -x]);
f=f0+pen*f;
