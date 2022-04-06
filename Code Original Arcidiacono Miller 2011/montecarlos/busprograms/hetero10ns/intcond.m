function [Like]=intcond(b,like,X);
global jjj

jjj=b;
U1=X*b;
p=exp(U1)./(1+exp(U1));
p=[p 1-p];

Like=-sum(log(max(sum(p.*like,2),eps)))./200;


