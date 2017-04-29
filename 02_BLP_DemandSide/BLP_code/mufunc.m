function f = mufunc(x2,theta2w)

global ns vfull 

[n k] = size(x2);
j = size(theta2w,2)-1;
mu = zeros(n,ns);

for i = 1:ns
      v_i = vfull(:,i:ns:k*ns);
 	  mu(:,i) = (x2.*v_i*theta2w(:,1));
end

f = mu;
