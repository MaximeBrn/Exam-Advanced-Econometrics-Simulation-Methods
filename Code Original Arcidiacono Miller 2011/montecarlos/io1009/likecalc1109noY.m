
function like=likecalc1109(b,Z2,FV,Firm2,N,Nf);



LikeD=((1-Firm2)+Firm2.*exp(Z2*b+FV))./(1+exp(Z2*b+FV));

%N refers to the number of rows in LikeY;

LikeD=reshape(LikeD,N,Nf);

like=prod(LikeD,2);
