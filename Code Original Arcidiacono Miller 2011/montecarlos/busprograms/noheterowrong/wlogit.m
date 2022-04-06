function [Like]=wlogit(b,Y,X,P);

U1=X*b;

Like=P'*(log(1+exp(U1))-Y.*U1);


