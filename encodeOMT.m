
function chromOMT=encodeOMT(E1)
kinds=size(E1,2);    
nums=size(E1,1); 
chromOMT=zeros(1,kinds);   %个体初始化

    for i=1:kinds
        chromOMT(i)=ceil(rand*nums); 
    end
end


