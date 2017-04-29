function [f]=gradient_vec(fcn,xx)

global mvalold mvalold0

h=10^(-6)*xx;

for i=1:size(xx,1)
    xxa=xx;
    xxb=xx;
    xxc=xx;
    xxd=xx;

    xxa(i,1)=xx(i,1)+h(i,1);
    xxb(i,1)=xx(i,1)-h(i,1);
    
    xxc(i,1)=xx(i,1)+2*h(i,1);
    xxd(i,1)=xx(i,1)-2*h(i,1);
    
    mvalold=mvalold0; f1f= feval(fcn,xxa); 
    mvalold=mvalold0; f1b= feval(fcn,xxb);    
    mvalold=mvalold0; f2f= feval(fcn,xxc);
    mvalold=mvalold0; f2b= feval(fcn,xxd);
    
    gradient_n(i,1)=(8*f1f-8*f1b-f2f+f2b)./(12*h(i,1));
end

f=gradient_n;
