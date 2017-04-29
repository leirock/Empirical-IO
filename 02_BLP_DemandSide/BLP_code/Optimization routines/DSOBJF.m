function f=dsobjf(x)
% Usage:
% f=dsobjf(x)
% Calculates the objective function value <f>
% at a point <x> for Shell Dual Problem
global A B C D E 
x=x(:); f=2*D'*x(1:5).^3 + (C*x(1:5))'*x(1:5) - B'*x(6:15);
