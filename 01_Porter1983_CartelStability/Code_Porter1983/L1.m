function f = L1(beta)
%This function computes the (-) maximum log likelihood
global n nb nb1d nb1s lambda se pn wght
load 'Data/jecnew';

%% Update lambda and wght series
if n==1; 
    load 'Data/wght0.csv'; %The 1st iteration weights come from a Probit estimation in Stata (see probit.do)
    wght=wght0;
    pn=wght;
    lambda=mean(pn);
else 
    lambda=lmda(beta);
end
    
%% Full-Information Maximum Likelihood Function (Davidson and MacKinnon, p524, eq.12.80)
GM=[1 -beta(nb1d+1+nb1s+1);-beta(nb1d+1) 1];
SG=[(beta(nb-2)) beta(nb);beta(nb) (beta(nb-1))]; %The Variance-Covariance matrix
B=[beta(1:nb1d)' 0 0 0 0 0; beta(nb1d+1+1) 0 beta(nb1d+1+2:nb1d+1+nb1s)']';
betadot=beta(1:nb1d+1+nb1s+1);
ydot=[Y(:,1);Y(:,2)];
W(:,size(W,2))=pn; %po replaced with the new weights
Z2(:,size(Z2,2))=pn; %po replaced with the new weights
Xdot=blkdiag([Z1 lngr],[Z2 lnQ]);
kost=-(G*T/2)*log(2*pi)-(T/2)*log(det(SG))+T*(log(abs(det(GM)))); %The first part of the negative log likelihood
f=-kost+.5*((ydot-Xdot*betadot)'*kron(inv(SG),eye(T))*(ydot-Xdot*betadot)); %The negative log likelihood

%% Calculate the standard errors
Ydot_new=W*B*(inv(GM)); %Y*GM=W*B+U -> Y=W*B*(GM^-1)+V
Xdot_new=blkdiag([Z1 Ydot_new(:,2)],[Z2 Ydot_new(:,1)]); %Xt=[Zt Yt], since y=X*beta+u=Z*beta1+Y*beta2+u
se=sqrt(diag(inv(Xdot_new'*(kron(inv(SG),eye(T)))*Xdot_new))); %Calculate the standard errors of estimaterd coefficients
se=[se;sqrt(diag(2.*kron(SG,SG)./T))]; %Add estimated standard errors to the vector se
se(length(se)-1)=[];  
%Estimated standard errors [shat1 shat12; shat12 shat2] -> [shat1 shat12 shat12 shat2]
%Delete one of the shat12 to make the vector dimension=37

end