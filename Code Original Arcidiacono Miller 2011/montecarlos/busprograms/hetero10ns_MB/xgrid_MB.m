
% Function to build transition matrix
function [xtran,xtranc]=xgrid(theta,xval);

% theta is a value of x2
% xval is all the possible values of x1t

n=length(xval); 
xub=[xval(2:n) ; inf]; % corresponds to x1_t+1
xtran=zeros(n,n);
xtranc=zeros(n,n);
lcdf=0;
for i=1:length(xval); % For each possible value of x1_t+1
        xtran(:,i)=(xub(i)>=xval).*(1-exp(-theta*(xub(i)-xval))-lcdf);
        lcdf=xtran(:,i)+lcdf;
        xtranc(:,i)=xtranc(:,i)+lcdf;
end;

