function [Like]=logitL(b,X);

U1=X*b;

Like=exp(U1)./(1+exp(U1));


