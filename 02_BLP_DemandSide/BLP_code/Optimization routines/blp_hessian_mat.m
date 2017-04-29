function hessian_mat = hessian(func,x,varargin)
global mvalold mvalold0

options_=[];
options_.gstep=0.001;
func = str2func(func);
n=size(x,1);
h1=max(abs(x),sqrt(options_.gstep)*ones(n,1))*eps^(1/6);
h_1=h1;
xh1=x+h1;
h1=xh1-x;
xh1=x-h_1;
h_1=x-xh1;
xh1=x;
  
mvalold=mvalold0;
f0=feval(func,x);
  
f1=zeros(size(f0,1),n);
f_1=f1;
for i=1:n    
    xh1(i)=x(i)+h1(i);

    mvalold=mvalold0;
    f1(:,i)=feval(func,xh1);

    xh1(i)=x(i)-h_1(i);

    mvalold=mvalold0;
    f_1(:,i)=feval(func,xh1);

    xh1(i)=x(i);
    i=i+1;
end
  
xh_1=xh1;
hessian_mat = zeros(size(f0,1),n*n);

for i=1:n    
    if i > 1
        k=[i:n:n*(i-1)];
        hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
    end

    hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
    temp=f1+f_1-f0*ones(1,n);
    for j=i+1:n
        xh1(i)=x(i)+h1(i);
        xh1(j)=x(j)+h_1(j);
        xh_1(i)=x(i)-h1(i);
        xh_1(j)=x(j)-h_1(j);

        mvalold=mvalold0;
        aa=feval(func,xh1);

        mvalold=mvalold0;
        bb=feval(func,xh_1);

        hessian_mat(:,(i-1)*n+j)=-(-aa-bb+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
        xh1(i)=x(i);
        xh1(j)=x(j);
        xh_1(i)=x(i);
        xh_1(j)=x(j);
        j=j+1;
    end
    i=i+1;
end 


