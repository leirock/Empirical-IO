function df = gradobj(theta2)
global invA1 IV1 ns cdid cdindex x2 v gmmresid mvalold

mval=mvalold;

theta2w = zeros(5,5);
theta2w(:,1)=theta2;

expmu = exp(mufunc(x2,theta2w));
shares = ind_sh(mval,expmu);
clear expmu

f1 = zeros(size(cdid,1),size(x2,2));

for ii=1:size(x2,2)
    xv = (x2(:,ii)*ones(1,ns)).*v(cdid,ns*(ii-1)+1:ns*ii);    
	temp = cumsum(xv.*shares);
	sum1 = temp(cdindex,:);
	sum1(2:size(sum1,1),:) = diff(sum1);
	f1(:,ii) = mean((shares.*(xv-sum1(cdid,:)))')';
	clear xv temp sum1
end
rel=size(theta2,1);

% computing (partial delta)/(partial theta2)
f = zeros(size(cdid,1),rel);
n = 1;
for i = 1:size(cdindex,1)
	temp = shares(n:cdindex(i),:);
	H1 = temp*temp';
	H = (diag(sum(temp')) - H1)/ns;
	f(n:cdindex(i),:) = - inv(H)*f1(n:cdindex(i),1:rel);
	n = cdindex(i) + 1;
end

df = 2*f'*IV1*invA1*IV1'*gmmresid;