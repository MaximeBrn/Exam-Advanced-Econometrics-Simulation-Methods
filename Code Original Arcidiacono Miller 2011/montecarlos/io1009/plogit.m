function [Like]=plogit(b,Y,X);

U1=X*b+adj;

prob=exp(U1)./(1+exp(U1));

Like=Y.*prob+(1-Y).*(1-prob);