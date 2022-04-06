function [Like,dg]=wlogit(b,Y,X,adj,P);

U1=X*b+adj;

Like=P'*(log(1+exp(U1))-Y.*U1);

if nargout>0

    dg=(P.*((1-(1./(1+exp(U1))))-Y))'*X;
    
end