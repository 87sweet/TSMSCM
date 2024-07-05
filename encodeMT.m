% function chromKMT=encodeKMT(C_KMT)
% d=size(C_KMT,2); 
% Z=zeros(NIND,d);
% Z(1,:)=rand(1,d);
% for i=2:NIND
%     Z(i,:)=4* Z(i-1,:).*(1-Z(i-1,:));
% end
% chromKMT=zeros(NIND,d);
% for i=1:NIND
%     chromKMT(i,:)=round(1+Z(i,:)*(d-1));
% end
function chromKMT=encodeKMT(C_KMT)
kinds=size(C_KMT,2);    
nums=size(C_KMT,1); 
chromKMT=zeros(1,kinds);   

    for i=1:kinds
        chromKMT(i)=ceil(rand*nums); 
%         sumW=sumW+C(chromKMT(i),i);
    end
end

