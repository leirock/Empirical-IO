function g=dscntg(x)
% Usage:
% g=dscntg(x)
% Calculates the gradient <g> of the constraint with the maximal
% residual at a point <x> for Shell Dual Problem
global A B C D E 
x=x(:); g=zeros(size(x));
[f,k]=max([A'*x(6:15) - 2*C*x(1:5) - 3*D.*x(1:5).^2 - E; -x]);
if f>0,  
    if k>5, g(k-5)=-1;
    else,   g(6:15)=A(:,k);
            g(1:5)=-2*C(:,k);
            g(k)=g(k)-6*D(k)*x(k);
    end 
end  
