
function [sumE_OMT,sumF_OMT,sumTime_OMT,sumCost_OMT]=chrom_EF_OMT(chromOMT,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back)
kinds=size(E1_OMT,2);    
sumE1_OMT=0;             
sumE2_OMT=0;             
sumE3_OMT=1;
sumE4_OMT=1;
sumE_OMT=0;
sumF_OMT=0;
sumTime_OMT=0;
sumCost_OMT=0;
u1=0.35;
u2=0.35;
u3=0.15;
u4=0.15;
for i=1:kinds
    sumTime_OMT=sumTime_OMT+T_back(chromOMT(i),i);
    sumCost_OMT=sumCost_OMT+C_cost_back(chromOMT(i),i);
    sumE1_OMT=sumE1_OMT+E1_OMT(chromOMT(i),i);
    sumE2_OMT=sumE2_OMT+E2_OMT(chromOMT(i),i);
    sumE3_OMT=sumE3_OMT*E3_OMT(chromOMT(i),i);
    sumE4_OMT=sumE4_OMT*E4_OMT(chromOMT(i),i);
    sumF_OMT=sumF_OMT+F_OMT(chromOMT(i),i);
%     fprintf('Individual %d: sumE_OMT = %.2f,sumE1_OMT = %.2f,sumE2_OMT = %.2f,sumE2_OMT = %.2f, sumE4_OMT = %.2f\n', i, sumE_OMT, sumE1_OMT,sumE2_OMT,sumE3_OMT,sumE4_OMT);
end
    sumE_OMT=sumE_OMT+u1*sumE1_OMT+u2*sumE2_OMT+u3*sumE3_OMT+u4*sumE4_OMT;
%     sumE_OMT=2*sumE_OMT;
end

