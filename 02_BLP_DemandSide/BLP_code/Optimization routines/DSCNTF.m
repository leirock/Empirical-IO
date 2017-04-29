function f=dscntf(x)
% Usage:
% f=dscntf(x)
% Calculates the maximal residual <f> for the set of constraints
% at a point <x> for Shell Dual Problem
global A B C D E
x=x(:);
f=max([A'*x(6:15) - 2*C*x(1:5) - 3*D.*x(1:5).^2 - E; -x]);
