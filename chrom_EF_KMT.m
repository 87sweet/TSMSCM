
function [sumE_KMT,sumF_KMT,sumTime_KMT,sumCost_KMT]=chrom_EF_KMT(chromKMT,E1_KMT,E2_KMT,E3_KMT,E4_KMT,F_KMT,T_front,C_cost_front)
kinds=size(E1_KMT,2);    
sumE1_KMT=0;              
sumE2_KMT=0;             
sumE3_KMT=1;
sumE4_KMT=1;
sumE_KMT=0;
sumF_KMT=0;
sumTime_KMT=0;
sumCost_KMT=0;
uT=0.35;
uC=0.35;
uRe=0.15;
uAv=0.15;
for i=1:kinds
    sumTime_KMT=sumTime_KMT+T_front(chromKMT(i),i);
    sumCost_KMT=sumCost_KMT+C_cost_front(chromKMT(i),i);
    sumE1_KMT=sumE1_KMT+E1_KMT(chromKMT(i),i);
    sumE2_KMT=sumE2_KMT+E2_KMT(chromKMT(i),i);
    sumE3_KMT=sumE3_KMT*E3_KMT(chromKMT(i),i);
    sumE4_KMT=sumE4_KMT*E4_KMT(chromKMT(i),i);
    sumF_KMT=sumF_KMT+F_KMT(chromKMT(i),i);
end
    sumE_KMT=sumE_KMT+uT*sumE1_KMT+uC*sumE2_KMT+uRe*sumE3_KMT+uAv*sumE4_KMT;
end

% function [sumE_KMT,sumF_KMT]=chrom_EF_KMT(chromKMT,E1_KMT,E2_KMT,E3_KMT,E4_KMT,F_KMT)
% kinds=size(E1_KMT,2);    %物品种类数目
% sumE_KMT=0;             %单个染色体的装包物品总重量 
% sumF_KMT=0;             %单个染色体的装包物品总重量
% uT=0.35;
% uC=0.35;
% uRe=0.15;
% uAv=0.15;
% for i=1:kinds
%     sumE_KMT=sumE_KMT+uT*E1_KMT(chromKMT(i),i)+uC*E2_KMT(chromKMT(i),i)+uRe*E3_KMT(chromKMT(i),i)+uAv*E4_KMT(chromKMT(i),i);
%     sumF_KMT=sumF_KMT+F_KMT(chromKMT(i),i);
% end
% end

