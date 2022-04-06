% addpath C:/Users/loicc/OneDrive/Documents/GitHub/Exam-Advanced-Econometrics-Simulation-Methods/Code Original Arcidiacono Miller 2011/montecarlos/io1009
cd 'C:/Users/loicc/OneDrive/Documents/GitHub/Exam-Advanced-Econometrics-Simulation-Methods/Code Original Arcidiacono Miller 2011/montecarlos/io1009'
clear

B=[];
T=[];
j=1;

while j<7
    
    if j==1

        load sim6f5sCom/bcom
    end
    if j==2
        load sim6f5sWrong/bwrong
    end
    if j==3
        load sim6f5sUpmodinit/bmodel
    end
    if j==4
        load sim6f5sUpdatainit/bdata
    end
    if j==5
        load sim6f5s2stage/bdata
    end
    if j==6
        load sim6f5snoY/bmodel
    end
    
    if j<3
        p=0;
    end
   
    
    if j>2
        p=exp(Binit(:,size(Binit,2)))./(1+exp(Binit(:,size(Binit,2))));
    end
    temp=[mean(Bccp)' std(Bccp)'];
    temp2=[mean(Byccp)' std(Byccp)';mean(p) std(p)];
    
    if j==2
        temp=[temp(1:2,:);0 0;temp(3:4,:)];
        temp2=[temp2(1:3,:);0 0; temp2(4:5,:)];
    end
    
    if j==6
        temp2=[zeros(5,2);mean(p) std(p)];
    end
    
    
    temp3=[temp;temp2];
    
    l=1;
    k=1;
    
    while k<12
        temp4(l:l+1,:)=[temp3(k,1);temp3(k,2)];
        k=k+1;
        l=l+2;
    end
    B=[B temp4];
    
    T=[T [mean(Tccp)./60;std(Tccp)./60]];
    j=j+1;
end

[B;T]