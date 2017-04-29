function f = var_cov(theta2)
global invA1 IV1 x1 x2 mvalold gmmresid

N = size(x1,1);
Z = size(IV1,2);

temp = jacob(mvalold,theta2);
a = [x1 temp]'*IV1;
IVres = IV1.*(gmmresid*ones(1,Z));

b = IVres'*IVres;

f = inv(a*invA1*a')*a*invA1*b*invA1*a'*inv(a*invA1*a');

