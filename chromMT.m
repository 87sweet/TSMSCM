

function [sumC_KMT,sumD_KMT]=chromCDKMT(chromKMT,C_KMT,D1_KMT,D2_KMT,D3_KMT)
kinds=size(C_KMT,2);    
sumC_KMT=0;              
sumD_KMT=0;            
uc=0.36;
ut=0.32;
ur=0.32;

for i=1:kinds
    sumC_KMT=sumC_KMT+C_KMT(chromKMT(i),i);
    sumD_KMT=sumD_KMT+uc*D1_KMT(chromKMT(i),i)+ut*D2_KMT(chromKMT(i),i)+ur*D3_KMT(chromKMT(i),i);
end
end

