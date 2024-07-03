
function chromKMT=encodeKMT(C_KMT)
kinds=size(C_KMT,2);    
nums=size(C_KMT,1); 
chromKMT=zeros(1,kinds); 

    for i=1:kinds
        chromKMT(i)=ceil(rand*nums); 
%         sumW=sumW+C(chromKMT(i),i);
    end
end

