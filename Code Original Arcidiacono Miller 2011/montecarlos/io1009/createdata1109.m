function [Firm,X,State,Y,Lfirm]=createdata1109(p,ctrans,S,Xn,Nf,Nm,T,bp,Tl)

Firm=zeros(Nm,T+Tl,Nf);
Lfirm=zeros(Nm,T+Tl+1,Nf);
X=unidrnd(Xn,Nm,1);
State=zeros(Nm,T+Tl+1);
Y=zeros(Nm,T+Tl);

State(:,1)=unidrnd(S,Nm,1);

Draw1=rand(Nm,T+Tl,Nf);
Draw2=rand(Nm,T+Tl);
Draw3=randn(Nm,T+Tl);


for nm=1:Nm
    Nfirm=0;
    for t=1:T+Tl
        
        for nf=1:Nf
            %given the number of others firms who were incumbents, the
            %state, and whether the firm was an incumbent, taking draws for
            %particular choice paths.
            Firm(nm,t,nf)=p(Nfirm-Lfirm(nm,t,nf)+1,X(nm),State(nm,t),Lfirm(nm,t,nf)+1)>Draw1(nm,t,nf);
            
        end
        
        Nfirm=sum(Firm(nm,t,:));
        Lfirm(nm,t+1,:)=Firm(nm,t,:);
        
        Y(nm,t)=[1 Nfirm X(nm)-1 State(nm,t)-1]*bp+Draw3(nm,t);
        
        State(nm,t+1)=1;
        
        for s=1:S-1
            
            State(nm,t+1)=State(nm,t+1)+(Draw2(nm,t)>ctrans(State(nm,t),s));
        end
        
    end
    
end

Firm=Firm(:,Tl+1:T+Tl,:);
State=State(:,Tl+1:T+Tl)-1;
Lfirm=Lfirm(:,Tl+1:T+Tl,:);
Y=Y(:,Tl+1:T+Tl);
