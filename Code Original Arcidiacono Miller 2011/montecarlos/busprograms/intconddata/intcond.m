function [Like]=intcond(b,like,X);

U1=X*b;
p=exp(U1)./(1+exp(U1));
p=[p 1-p];

Like=-sum(log(sum(p.*like,2)));


