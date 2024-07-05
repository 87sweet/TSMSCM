% function chromKMT=encodeOMT(C_KMT)
% d=size(E1_OMT,2); 
% Z=zeros(NINDO,d);
% Z(1,:)=rand(1,d);
% for i=2:NIND
%     Z(i,:)=4* Z(i-1,:).*(1-Z(i-1,:));
% end
% chromOMT=zeros(NIND,d);
% for i=1:NIND
%     chromOMT(i,:)=round(1+Z(i,:)*(d-1));
% end
function chromOMT=encodeOMT(E1)
kinds=size(E1,2);    
nums=size(E1,1); 
chromOMT=zeros(1,kinds);   %个体初始化

    for i=1:kinds
        chromOMT(i)=ceil(rand*nums); 
    end
end


