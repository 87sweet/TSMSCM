

function [sumC_OMT,sumD_OMT]=chrom_CD_OMT(chrom0MT,C_OMT,D1_OMT,D2_OMT,D3_OMT)
kinds=size(C_OMT,2);   
sumC_OMT=0;              
sumD_OMT=0;             
uc=0.36;
ut=0.32;
ur=0.32;

for i=1:kinds
    sumC_OMT=sumC_OMT+C_OMT(chrom0MT(i),i);
    sumD_OMT=sumD_OMT+uc*D1_OMT(chromOMT(i),i)+ut*D2_OMT(chromOMT(i),i)+ur*D3_OMT(chromOMT(i),i);
end
end

