

function ChromOMT=InitPopOMT(NINDO,E1)
N=size(E1,2);    %��Ʒ������Ŀ
ChromOMT=zeros(NINDO,N);%���ڴ洢��Ⱥ
for i=1:NINDO
    ChromOMT(i,:)=encodeKMT(E1);%������ɳ�ʼ��Ⱥ
end
end