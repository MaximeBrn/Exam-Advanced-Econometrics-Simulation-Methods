

function [Fv,BigP]=match1109(NFirm2,LFirm2,X2,State2,fv,bigp,Nobs,N);

Fv=zeros(Nobs,1);
BigP=zeros(Nobs,N);
n=1;

while n<Nobs+1
    
    Fv(n,1)=fv(NFirm2(n)+1,X2(n),State2(n),LFirm2(n)+1);
    
    BigP(n,:)=bigp(NFirm2(n)+1,:,X2(n),State2(n),LFirm2(n)+1);
    
    n=n+1;
    
end