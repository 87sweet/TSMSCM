
function ChromKMT=InitPopKMT(NIND,C_KMT)
N=size(C_KMT,2);    %物品种类数目
ChromKMT=zeros(NIND,N);%用于存储种群
for i=1:NIND
    ChromKMT(i,:)=encodeKMT(C_KMT);%随机生成初始种群
end
end