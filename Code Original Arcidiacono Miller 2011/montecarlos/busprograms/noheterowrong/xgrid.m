function [xtran,xtranc]=xgrid(theta,xval);

n=length(xval);
xub=[xval(2:n) ; inf];
xtran=zeros(n,n);
xtranc=zeros(n,n);
lcdf=0;
for i=1:length(xval);
        xtran(:,i)=(xub(i)>=xval).*(1-exp(-theta*(xub(i)-xval))-lcdf);
        lcdf=xtran(:,i)+lcdf;
        xtranc(:,i)=xtranc(:,i)+lcdf;
end;

