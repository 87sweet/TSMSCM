
function allObj_EF_KMT=allObject_EF_KMT(chromKMT,E1_KMT,E2_KMT,E3_KMT,E4_KMT,F_KMT,T_front,C_cost_front)
b=size(chromKMT,1); %种群大小
allObj_EF_KMT=zeros(b,5);
w21=0.67;
w22=0.33;
for i=1:b
        [sumE_KMT,sumF_KMT,sumTime_KMT,sumCost_KMT]=chrom_EF_KMT(chromKMT(i,:),E1_KMT,E2_KMT,E3_KMT,E4_KMT,F_KMT,T_front,C_cost_front);
        allObj_EF_KMT(i,1)=sumTime_KMT;
        allObj_EF_KMT(i,2)=sumCost_KMT;
        allObj_EF_KMT(i,3)=sumE_KMT;
        allObj_EF_KMT(i,4)=sumF_KMT;
        allObj_EF_KMT(i,5)=w21*sumE_KMT+w22*sumF_KMT;
end 
end



