function f = LL0(beta)
%This function computes the (-) maximum log likelihood

global nb nb1d nb1s
load 'Data/jecnew0';

GM=[1 -beta(nb1d+1+nb1s);-beta(nb1d+1) 1];
SG=[(beta(nb-2)) beta(nb);beta(nb) (beta(nb-1))];
B=[beta(1:nb1d)' 0 0 0 0; beta(nb1d+1+1) 0 beta(nb1d+1+2:nb1d+1+nb1s)']';
betadot=beta(1:nb1d+1+nb1s+1);
ydot=[Y(:,1);Y(:,2)];
Xdot=blkdiag([Z1 lngr],[Z2 lnQ]);

kost=-(G*T/2)*log(2*pi)-(T/2)*log(det(SG))+T*(log(abs(det(GM)))); %The first part of the log likelihood
f=-kost+.5*((ydot-Xdot*betadot)'*kron(inv(SG),eye(T))*(ydot-Xdot*betadot)); %(-) the log likelihood