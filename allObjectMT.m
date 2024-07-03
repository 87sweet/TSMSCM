
function allObjMT=allObjectMT(ChromKMT,ChromOMT,C,D1,D2,D3,E1,E2,E3,E4,F)
Chromtotal = [ChromKMT, ChromOMT];
a=size(Chromtotal,1); %种群大小
allObjtotal=zeros(a,1);
w1=0.58;
w2=0.42;
kinds=size(C_KMT,2);    %物品种类数目
sumC_KMT=0;             %单个染色体的装包物品总重量 
sumD_KMT=0;             %单个染色体的装包物品总重量
uc=0.36;
ut=0.32;
ur=0.32;

for i=1:kinds
    sumC_KMT=sumC_KMT+C_KMT(chromKMT(i),i);
    sumD_KMT=sumD_KMT+uc*D1_KMT(chromKMT(i),i)+ut*D2_KMT(chromKMT(i),i)+ur*D3_KMT(chromKMT(i),i);
end

for i=1:a
        %[sumC_KMT,sumD_KMT]=chromCDKMT(ChromKMT(i,:),C_KMT,D_KMT);
        [sumC_KMT,sumD_KMT]=chromCDKMT(ChromKMT(i,:),C_KMT,D1_KMT,D2_KMT,D3_KMT);
        allObjKMT(i,:)=w1*sumC_KMT+w2*sumD_KMT;
end 
end



