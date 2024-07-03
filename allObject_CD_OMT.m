
function allObj_CD_OMT=allObject_CD_OMT(ChromOMT,C_OMT,D1_OMT,D2_OMT,D3_OMT)
a=size(ChromOMT,1); %种群大小
allObj_CD_OMT=zeros(a,3);
w1=0.58;
w2=0.42;
for i=1:a
        %[sumC_KMT,sumD_KMT]=chromCDKMT(ChromKMT(i,:),C_KMT,D_KMT);
        [sumC_OMT,sumD_OMT]=chrom_CD_OMT(ChromOMT(i,:),C_OMT,D1_OMT,D2_OMT,D3_OMT);
        allObj_CD_OMT(i,1)=sumC_OMT;
        allObj_CD_OMT(i,2)=sumD_OMT;
        allObj_CD_OMT(i,3)=w1*sumC_OMT+w2*sumD_OMT;
end 
end



