function [f,g] = gmmobj3(theta2,par)
global theta1 x1 IV1 invA1 gmmresid ppp mymaxfunevals...
fcnevals fvals_track

par=0;

if size(theta2,1)==1
    theta2=theta2';
end

delta = meanval(theta2);

if max(isnan(delta)) == 1
	f=1e+10;	
    gmmresid=(1e+10)*ones(size(delta));
    g=(1e+10)*ones(size(theta2));         
    fcnevals=fcnevals+1;
    if fcnevals<=mymaxfunevals
        fvals_track(fcnevals,ppp)=f;
    end        
else
	temp1 = x1'*IV1;
	temp2 = delta'*IV1;
    theta1 = inv(temp1*invA1*temp1')*temp1*invA1*temp2';
	gmmresid = delta - x1*theta1;
	temp1 = gmmresid'*IV1;
	f = temp1*invA1*temp1';
    g=gradobj(theta2);
    fcnevals=fcnevals+1;
    if fcnevals<=mymaxfunevals
        fvals_track(fcnevals,ppp)=f;
    end    
end
