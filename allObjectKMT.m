

function allObjKMT=allObjectKMT(ChromKMT,C,D)
a=size(ChromKMT,1); %种群大小
allObjKMT=zeros(a,1);
w1=1;
w2=1;
for i=1:a
        [sumC,sumD]=chromCDKMT(ChromKMT(i,:),C,D);
        allObjKMT(i,:)=w1*sumC+w2*sumD;
end 
end



