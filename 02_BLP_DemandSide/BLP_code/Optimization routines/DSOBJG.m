function g=dsobjg(x)
% Usage:
% g=dsobjg(x)
% Calculates the gradient <g> of the objective function
% at a point <x> for Shell Dual Problem
global A B C D E 
x=x(:);g=x; 
g(1:5)=6*D.*x(1:5).^2 + 2*C*x(1:5);
g(6:15)= - B;
