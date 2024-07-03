
function allObj_CD_KMT=allObject_CD_KMT(ChromKMT,C_KMT,D1_KMT,D2_KMT,D3_KMT)
a=size(ChromKMT,1); %种群大小
allObj_CD_KMT=zeros(a,3);
w1=0.58;
w2=0.42;
for i=1:a
        %[sumC_KMT,sumD_KMT]=chromCDKMT(ChromKMT(i,:),C_KMT,D_KMT);
        [sumC_KMT,sumD_KMT]=chrom_CD_KMT(ChromKMT(i,:),C_KMT,D1_KMT,D2_KMT,D3_KMT);
        allObj_CD_KMT(i,1)=sumC_KMT;
        allObj_CD_KMT(i,2)=sumD_KMT;
        allObj_CD_KMT(i,3)=w1*sumC_KMT+w2*sumD_KMT;
end 
end



