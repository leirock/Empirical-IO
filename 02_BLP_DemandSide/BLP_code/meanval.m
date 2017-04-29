function f = meanval(theta2)

global x2 s_jt mvalold oldt2 mvalolds mval_track

if max(abs(theta2-oldt2)) < 0.01;
	tol = 1e-9;
	flag = 0;
else
  	tol = 1e-6;
	flag = 1;
end

theta2w = zeros(5,5); theta2w(:,1)=theta2;
expmu = exp(mufunc(x2,theta2w));
norm = 1;
avgnorm = 1;

i = 0;
mvalolds=[mvalold];

while norm > tol*10^(flag*floor(i/50)) && avgnorm > 1e-3*tol*10^(flag*floor(i/50))
%while (norm > 1e-16 && i<=5000)
	  mval = mvalold.*s_jt./mktsh(mvalold,expmu); 
      t = abs(mval-mvalold);
	  norm = max(t);
      avgnorm = mean(t);
  	  mvalold = mval;
      i = i + 1;
end

if flag == 1 && max(isnan(mval)) < 1;
   mvalold = mval;
   oldt2 = theta2;
end

f = log(mval);

mvalolds=[mvalolds,mvalold];
mval_track=[mval_track,mvalold];
