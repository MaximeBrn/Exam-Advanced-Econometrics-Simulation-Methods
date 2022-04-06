function [Like]=logitL(b,Y,X);

U1=X*b;

Like=(Y.*exp(U1)+(1-Y))./(1+exp(U1));


