
function like=likecalc1109(bp,sig,b,Y,ZZ2,Z2,FV,Firm2,N,Nf);



EY=ZZ2*bp;

LikeY=normpdf(Y-EY,0,sig);

LikeD=((1-Firm2)+Firm2.*exp(Z2*b+FV))./(1+exp(Z2*b+FV));

%N refers to the number of rows in LikeY;

LikeD=reshape(LikeD,N,Nf);

like=prod(LikeD,2).*LikeY;
