

function ChromOMT=InitPopOMT(NINDO,E1)
N=size(E1,2);    %物品种类数目
ChromOMT=zeros(NINDO,N);%用于存储种群
for i=1:NINDO
    ChromOMT(i,:)=encodeKMT(E1);%随机生成初始种群
end
end