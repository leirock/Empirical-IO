function f = mktsh(mval, expmu)

global ns 

f = sum((ind_sh(mval,expmu))')/ns;

f = f';