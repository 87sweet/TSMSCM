
function ChromKMT=InitPopKMT(NIND,C_KMT)
N=size(C_KMT,2);    %��Ʒ������Ŀ
ChromKMT=zeros(NIND,N);%���ڴ洢��Ⱥ
for i=1:NIND
    ChromKMT(i,:)=encodeKMT(C_KMT);%������ɳ�ʼ��Ⱥ
end
end