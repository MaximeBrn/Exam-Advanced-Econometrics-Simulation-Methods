function [PType]=intcondP(b,like,X);

U1=X*b;
p=exp(U1)./(1+exp(U1));
p=[p 1-p];

PType=(like.*p)./(sum(like.*p,2)*ones(1,2));




