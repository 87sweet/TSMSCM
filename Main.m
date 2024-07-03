clear;clc;format compact;tic;hold on

% �������ݷ�Χ
Co_range = [0.85, 0.95];
D_range = [0.85, 0.95];
T_range = [15, 25];
C_cost_range1 = [2000, 3000];
C_cost_range2 = [1000, 1500];
% C_cost_range1 = [1000, 1500];
% C_cost_range2 = [1000, 1500];
Re_range = [0.9, 1];
Av_range = [0.9, 1];
Sta_range = [0.85, 0.95];

% ����ǰ3�кͺ�9�е�����
front_rows = 30;
back_rows = 30;

aa=3;
bb=9;

% ���ɷ�����̬�ֲ����������
rng('default'); % For reproducibility

Co_front = normrnd(mean(Co_range), 0.05, [front_rows, aa]);
Co_back = normrnd(mean(Co_range), 0.05, [back_rows, bb]);

D1_front = normrnd(mean(D_range), 0.05, [front_rows, aa]);
D1_back = normrnd(mean(D_range), 0.05, [back_rows, bb]);
D2_front = normrnd(mean(D_range), 0.05, [front_rows, aa]);
D2_back = normrnd(mean(D_range), 0.05, [back_rows, bb]);
D3_front = normrnd(mean(D_range), 0.05, [front_rows, aa]);
D3_back = normrnd(mean(D_range), 0.05, [back_rows, bb]);

T_front = randi(T_range, [front_rows, aa]);
T_back = randi(T_range, [back_rows, bb]);

C_cost_front = randi(C_cost_range1, [front_rows, aa]);
C_cost_back = randi(C_cost_range2, [back_rows, bb]);

Re_front = normrnd(mean(Re_range), 0.02, [front_rows, aa]);
Re_back = normrnd(mean(Re_range), 0.02, [back_rows, bb]);

Av_front = normrnd(mean(Av_range), 0.02, [front_rows, aa]);
Av_back = normrnd(mean(Av_range), 0.02, [back_rows, bb]);

Sta_front = normrnd(mean(Sta_range), 0.05, [front_rows, aa]);
Sta_back = normrnd(mean(Sta_range), 0.05, [back_rows, bb]);

% �ϲ�ǰ3�кͺ�9�е�����
Co = [Co_front, Co_back];
D1 = [D1_front, D1_back];  %F1
D2 = [D2_front, D2_back];  %F2
D3 = [D3_front, D3_back];  %F3
T = [T_front, T_back];  %T
C = [C_cost_front, C_cost_back];  %C
Re = [Re_front, Re_back];  %Re
Av = [Av_front, Av_back];  %Av
Sta = [Sta_front, Sta_back];  %Sta

% ���������С��һ������
%min_max_normalize = @(x) (x - min(x(:))) / (max(x(:)) - min(x(:)));
max_min_normalize = @(x) (max(x(:))-x) / (max(x(:)) - min(x(:)));

T_normalized = max_min_normalize(T);
C_normalized = max_min_normalize(C);

% ��������
C_KMT = Co(:, 1:aa);
C_OMT = Co(:, aa+1:end);

D1_KMT = D1(:, 1:aa);
D2_KMT = D2(:, 1:aa);
D3_KMT = D3(:, 1:aa);
D1_OMT = D1(:, aa+1:end);
D2_OMT = D2(:, aa+1:end);
D3_OMT = D3(:, aa+1:end);

E1=T_normalized;
E2=C_normalized;
E3=Re;
E4=Av;
E1_KMT = T_normalized(:, 1:aa);
E2_KMT = C_normalized(:, 1:aa);
E3_KMT = Re(:, 1:aa);
E4_KMT = Av(:, 1:aa);
E1_OMT = T_normalized(:, aa+1:end);
E2_OMT = C_normalized(:, aa+1:end);
E3_OMT = Re(:, aa+1:end);
E4_OMT = Av(:, aa+1:end);
F = Sta;
F_KMT = Sta(:, 1:aa);
F_OMT = Sta(:, aa+1:end);

% ����Ϊ.mat�ļ�
% save('data.mat', 'C_KMT', 'C_OMT', 'D1_KMT', 'D2_KMT', 'D3_KMT', 'D1_OMT', 'D2_OMT', 'D3_OMT', 'E1_KMT', 'E2_KMT', 'E3_KMT', 'E4_KMT', 'E1', 'E2', 'E3', 'E4', 'F_KMT', 'F');
save('data.mat', 'Co_range', 'D_range', 'T_range', 'C_cost_range1', 'C_cost_range2', ...
    'Re_range', 'Av_range', 'Sta_range', 'front_rows', 'back_rows', 'aa', 'bb', ...
    'Co_front', 'Co_back', 'D1_front', 'D1_back', 'D2_front', 'D2_back', ...
    'D3_front', 'D3_back', 'T_front', 'T_back', 'C_cost_front', 'C_cost_back', ...
    'Re_front', 'Re_back', 'Av_front', 'Av_back', 'Sta_front', 'Sta_back', ...
    'Co', 'D1', 'D2', 'D3', 'T', 'C', 'Re', 'Av', 'Sta', ...
    'T_normalized', 'C_normalized',...
    'C_KMT','C_OMT','D1_KMT','D2_KMT', 'D3_KMT', 'D1_OMT', 'D2_OMT', 'D3_OMT', 'E1_KMT', 'E2_KMT', 'E3_KMT', 'E4_KMT', 'E1_OMT', 'E2_OMT', 'E3_OMT', 'E4_OMT', 'F_KMT', 'F_OMT');

load data.mat
%cap=100;

%% ��������
NIND=20;        %��Ⱥ��С
NINDO=20;        %��Ⱥ��С 
MAXGEN=100;     %��������
Pc=0.9;         %�������
Pm=0.09;        %�������
%N=size(C,2);

%% ��ʼ����Ⱥ 

%% �Ż�
N=1;
ChromKMT_xnew1=InitPopKMT(NIND,C_KMT);
xnew1=ChromKMT_xnew1;
allObj_CD_KMT_xnew1=allObject_CD_KMT(ChromKMT_xnew1,C_KMT,D1_KMT,D2_KMT,D3_KMT); 
[value,index]=max(allObj_CD_KMT_xnew1(:,3)); 
tourPbest=ChromKMT_xnew1; 
tourGbest=ChromKMT_xnew1(index,:); 
recordPbest=-inf*ones(NIND,1);               %�������ż�¼
recordGbest=allObj_CD_KMT_xnew1(index,3);                %Ⱥ�����ż�¼
L_best=zeros(MAXGEN,1);          %��¼ÿ�ε���������ȫ�����Ÿ�����ܾ���  %L_best(N)=allObj_CD_KMT(index1);
L_average=zeros(MAXGEN,size(allObj_CD_KMT_xnew1, 2));

ChromKMT=InitPopKMT(NIND,C_KMT);
allObj_CD_KMT=allObject_CD_KMT(ChromKMT,C_KMT,D1_KMT,D2_KMT,D3_KMT);
for N=1:MAXGEN  
    ChromOMT_xnew2=InitPopOMT(NINDO,E1_OMT);
    xnew2=ChromOMT_xnew2;
    allObj_EF_OMT_xnew2=allObject_EF_OMT(ChromOMT_xnew2,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back); 
    [value111,index111]=max(allObj_EF_OMT_xnew2(:,5)); 
    tourPbest0=ChromOMT_xnew2; 
    tourGbestO=ChromOMT_xnew2(index111,:); 
    recordPbestO=-inf*ones(NINDO,1);               %�������ż�¼
    recordGbestO=allObj_EF_OMT_xnew2(index111,5);                %Ⱥ�����ż�¼
    L_bestO=zeros(MAXGEN,1);          %��¼ÿ�ε���������ȫ�����Ÿ�����ܾ���  %L_best(N)=allObj_CD_KMT(index1);
    L_averageO=zeros(MAXGEN,size(allObj_EF_OMT_xnew2, 2));
    
    ChromOMT=InitPopOMT(NINDO,E1_OMT);
    allObj_EF_OMT=allObject_EF_OMT(ChromOMT,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back); 
    

    for M=1:MAXGEN
        for i=1:NINDO
            allObj_EF_KMT=allObject_EF_KMT(ChromKMT,E1_KMT,E2_KMT,E3_KMT,E4_KMT,F_KMT,T_front,C_cost_front);
            
            a = 250; % ʾ��ֵ�������ʵ���������/ʱ��Լ��
            b = 16000; % ʾ��ֵ�������ʵ���������/�ɱ�Լ��

            % ��ȡ allObj_EF_OMT �� allObj_EF_KMT �����е�ֵ
            OMT_col_time = allObj_EF_OMT(:, 1);
            KMT_col_time = allObj_EF_KMT(:, 1);
            total_time=OMT_col_time+KMT_col_time;
            OMT_col_cost = allObj_EF_OMT(:, 2);
            KMT_col_cost = allObj_EF_KMT(:, 2);
            total_cost=OMT_col_cost+ KMT_col_cost;
            % ���ÿ�е�ֵ�Ƿ�����Լ������
             for i = 1:length(OMT_col_time)
                 if total_time(i) > a
                     allObj_EF_OMT(i, 5) = 0;
                 end
                 if total_cost(i) > b
                 allObj_EF_OMT(i, 5) = 0;
                 end
             end

             if allObj_EF_OMT(i,5) > recordPbestO(i)
                 recordPbestO(i) = allObj_EF_OMT(i,5);
                 tourPbestO(i,:) = ChromOMT(i,:);
             end
             if allObj_EF_OMT(i,5) > recordGbestO
                 recordGbestO = allObj_EF_OMT(i,5);
                 tourGbestO = ChromOMT(i,:);
             end
                       
             % 1111��������Ž��н���
            L=size(E1_OMT,2);
            c1=randi([1 L],1,1); %��������λ
            c2=randi([1 L],1,1); %��������λ
            while c1==c2
                c2=randi([1 L],1,1);
            end
            
                        % ����͸��Ƴ�ʼ����Ԫ��
            a0 = xnew2(i,:);  % ����xnew2�ĵ�i�е�a0���������ݺͲο�
            b0 = tourPbestO(i,:);  % ����tourPbestO�ĵ�i�е�b0

            % ȷ�������ı߽�
            chb1 = min(c1, c2);  % ȷ���������½�
            chb2 = max(c1, c2);  % ȷ���������Ͻ�

            % ����ѭ��
            for j = chb1:chb2
                % ��xnew2�ĵ�i�е�j�з���b0�ĵ�j��Ԫ��
                xnew2(i,j) = b0(j);  
            end

            % ���Ҳ��������ܵ��ظ�Ԫ��
            for j = chb1:chb2
                x = find(xnew2(i,:) == xnew2(i,j));  % ����xnew2��i�����뵱ǰ����λ��Ԫ����ȵ�����λ��
                i1 = x(x ~= j);  % ���ҵ���λ�����ų���ǰ������λ��
                if ~isempty(i1)  % ��������ظ���Ԫ��
                    % ��ԭʼ��a0�е�Ԫ���滻���ظ�Ԫ�أ�ȷ��ÿ��Ԫ����Ψһ��
                    xnew2(i,i1) = a0(j);  
                end
            end
            
            dist=allObject_EF_OMT(xnew2,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back); %allObjectKMT ��Ŀ��ĺ���
            if  allObj_EF_OMT(i,5)<dist(i,5) %allObjKMTΪĿ�꺯��allObjectKMT�󵽵�Ŀ��ֵ
            ChromOMT(i,:)=xnew2(i,:);
            end 
            
            % �������
            c1 = randi([1 L], 1, 1); % ��������λ
            c2 = randi([1 L], 1, 1); % ��������λ

            if c1 ~= c2
                a0 = xnew2(i, :);  % ����ԭʼ�����е�״̬
                b0 = tourGbestO;   % ��ȡȺ�����Ž�

                chb1 = min(c1, c2);  % ȷ���������½�
                chb2 = max(c1, c2);  % ȷ���������Ͻ�

                % ���н������
                for j = chb1:chb2
                    xnew2(i, j) = b0(j);  % ��b0�ĵ�j��Ԫ�ظ��Ƶ�xnew2�ĵ�i�е�j��
                end

                % ���Ҳ��������ܵ��ظ�Ԫ��
                for j = chb1:chb2
                    x = find(xnew2(i, :) == xnew2(i, j));  % ����xnew2��i�����뵱ǰ����λ��Ԫ����ȵ�����λ��
                    i1 = x(x ~= j);  % ���ҵ���λ�����ų���ǰ������λ��
                    if ~isempty(i1)  % ��������ظ���Ԫ��
                        xnew2(i, i1) = a0(j);  % ��ԭʼ��a0�е�Ԫ���滻���ظ�Ԫ��
                    end
                end
            end

            % �����½��Ŀ�꺯��ֵ���Ƚ�
            dist2 = allObject_EF_OMT(xnew2, E1_OMT, E2_OMT, E3_OMT, E4_OMT, F_OMT, T_back, C_cost_back);
            if allObj_EF_OMT(i) < dist2
                ChromOMT(i, :) = xnew2(i, :);
            end
            
        
        %% �������
            c1=randi([1 L],1,1); %��������λ
            c2=randi([1 L],1,1); %��������λ
            while c1==c2
                c2=randi([1 L],1,1);
            end
            temp=xnew2(i,c1);
            xnew2(i,c1)=xnew2(i,c2);
            xnew2(i,c2)=temp;
        
            dist3=allObject_EF_OMT(xnew2,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back);
            if  allObj_EF_OMT(i)<dist3
            ChromOMT(i,:)=xnew2(i,:);
            end
        end
        
        for i=1:NINDO
            allObj_EF_OMT=allObject_EF_OMT(ChromOMT,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back);  
            allObj_EF_KMT=allObject_EF_KMT(ChromKMT,E1_KMT,E2_KMT,E3_KMT,E4_KMT,F_KMT,T_front,C_cost_front);
            % ��ȡ allObj_EF_OMT �� allObj_EF_KMT �����е�ֵ
            OMT_col_time = allObj_EF_OMT(:, 1);
            KMT_col_time = allObj_EF_KMT(:, 1);
            total_time=allObj_EF_OMT(:, 1)+ allObj_EF_KMT(:, 1);
            OMT_col_cost = allObj_EF_OMT(:, 2);
            KMT_col_cost = allObj_EF_KMT(:, 2);
            total_cost=allObj_EF_OMT(:, 2)+ allObj_EF_KMT(:, 2);
        end
        
        for i = 1:NINDO
            if allObj_EF_OMT(i,5) > recordPbestO(i)
                 recordPbestO(i) = allObj_EF_OMT(i,5);
                 tourPbestO(i,:) = ChromOMT(i,:);
            end
            if allObj_EF_OMT(i,5) > recordGbestO
                  recordGbestO = allObj_EF_OMT(i,5);
                  tourGbestO = ChromOMT(i,:);
            end
        end 
        L_bestO(M) = recordGbestO;
        L_averageO(M,:) = mean(allObj_EF_OMT, 1);
        fprintf('Iteration %d: recordGbestO = %f, tourGbestO = %s\n', i, recordGbestO, mat2str(tourGbestO));
    end

    
%% ���²�ģ�����Ž������ϲ�ģ�ͣ�����ϲ�ģ�͵����Ž�
    for i=1:NIND
        allObj_CD_KMT=allObject_CD_KMT(ChromKMT,C_KMT,D1_KMT,D2_KMT,D3_KMT);  
        allObj_CD_OMT=allObject_CD_KMT(tourGbestO,C_OMT,D1_OMT,D2_OMT,D3_OMT);       
        c = 10.6; % ʾ��ֵ�������ʵ���������/ʱ��Լ��

        % ��ȡ allObj_EF_OMT �� allObj_EF_KMT �����е�ֵ
        OMT_col_Co = allObj_CD_OMT(:,1);
        KMT_col_Co = allObj_CD_KMT(:,1);
        total_Co=allObj_CD_KMT(:, 1)+allObj_CD_OMT(:,1);
        % ���ÿ�е�ֵ�Ƿ���������
        for i = 1:length(OMT_col_Co)
            if total_Co(i)< c
                allObj_CD_OMT(i, 3) = 0;
            end
        end
        
        % 1111��������Ž��н���
        L=size(C_KMT,2);
        c1=randi([1 L],1,1); %��������λ
        c2=randi([1 L],1,1); %��������λ
        while c1==c2
            c2=randi([1 L],1,1);
        end
            
         % ����͸��Ƴ�ʼ����Ԫ��
         a0 = xnew1(i,:);  % ����xnew2�ĵ�i�е�a0���������ݺͲο�
         b0 = tourPbest(i,:);  % ����tourPbestO�ĵ�i�е�b0

         % ȷ�������ı߽�
         chb1 = min(c1, c2);  % ȷ���������½�
         chb2 = max(c1, c2);  % ȷ���������Ͻ�

         % ����ѭ��
         for j = chb1:chb2
             % ��xnew2�ĵ�i�е�j�з���b0�ĵ�j��Ԫ��
              xnew1(i,j) = b0(j);  
         end

            % ���Ҳ��������ܵ��ظ�Ԫ��
         for j = chb1:chb2
              x = find(xnew1(i,:) == xnew1(i,j));  % ����xnew2��i�����뵱ǰ����λ��Ԫ����ȵ�����λ��
              i1 = x(x ~= j);  % ���ҵ���λ�����ų���ǰ������λ��
              if ~isempty(i1)  % ��������ظ���Ԫ��
                    % ��ԭʼ��a0�е�Ԫ���滻���ظ�Ԫ�أ�ȷ��ÿ��Ԫ����Ψһ��
                  xnew1(i,i1) = a0(j);  
              end
          end
            
          dist=allObject_CD_KMT(xnew1,C_KMT,D1_KMT,D2_KMT,D3_KMT); %allObjectKMT ��Ŀ��ĺ���
          if  allObj_CD_KMT(i,3)<dist(i,3) %allObjKMTΪĿ�꺯��allObjectKMT�󵽵�Ŀ��ֵ
          ChromKMT(i,:)=xnew1(i,:);
          end 

       % 2222��Ⱥ�����Ž��н���
        c1=randi([1 L],1,1); %��������λ
        c2=randi([1 L],1,1); %��������λ
        if c1~=c2
            a0=xnew1(i,:);b0=tourGbest; 
            chb1=min(c1,c2);     
            chb2=max(c1,c2);
            for j=chb1:chb2
                a1=xnew1(i,:); b1=tourGbest;
                xnew1(i,j)=b0(j);  
                x=find(xnew1(i,:)==xnew1(i,j));  %���ҵ�i������a(j)ֵ��ȵ�Ԫ�ص�λ�ã������������������x�С�
                i1=x(x~=j); %������x���ų�λ��Ϊj��Ԫ�أ����������������i1�С�
                if ~isempty(i1)  %�ж�����i1�Ƿ�Ϊ�ա�
                xnew1(i,i1)=a1(j); %�������i1�ǿգ�������a��λ��Ϊi1��Ԫ�ظ���Ϊ����a1�еĵ�i��Ԫ�ء�
                end
            end
        end
        dist2=0;     
        dist2=dist2+allObject_CD_KMT(xnew1,C_KMT,D1_KMT,D2_KMT,D3_KMT);
        if  allObj_CD_KMT<dist2
            ChromKMT(i,:)=xnew1(i,:);
        end
        
        %% �������
        c1=randi([1 L],1,1); %��������λ
        c2=randi([1 L],1,1); %��������λ
        while c1==c2
%             c1=randi([1 L],1,1); %��������λ
            c2=randi([1 L],1,1);
        end
        temp=xnew1(i,c1);
        xnew1(i,c1)=xnew1(i,c2);
        xnew1(i,c2)=temp;
        
        dist3=allObject_CD_KMT(xnew1,C_KMT,D1_KMT,D2_KMT,D3_KMT);
        if  allObj_CD_KMT(i,3)<dist3(i,3)
           ChromKMT(i,:)=xnew1(i,:);
        end
    end
    
    for i=1:NIND
        allObj_CD_KMT=allObject_CD_KMT(ChromKMT,C_KMT,D1_KMT,D2_KMT,D3_KMT);  
        allObj_CD_OMT=allObject_CD_KMT(tourGbestO,C_OMT,D1_OMT,D2_OMT,D3_OMT);       
        c = 11; % ʾ��ֵ�������ʵ���������/ʱ��Լ��

        % ��ȡ allObj_EF_OMT �� allObj_EF_KMT �����е�ֵ
        OMT_col_Co = allObj_CD_OMT(:,1);
        KMT_col_Co = allObj_CD_KMT(:,1);
        total_Co=allObj_CD_KMT(:, 1)+allObj_CD_OMT(:,1);
        % ���ÿ�е�ֵ�Ƿ���������
        for i = 1:length(OMT_col_Co)
            if total_Co(i)< c
                allObj_CD_KMT(i, 3) = 0;
            end
        end
    end
   
    allObj_CD_KMT = allObject_CD_KMT(ChromKMT, C_KMT, D1_KMT, D2_KMT, D3_KMT);
    for i = 1:NIND
        if allObj_CD_KMT(i,3) > recordPbest(i)
            recordPbest(i) = allObj_CD_KMT(i,3);
            tourPbest(i,:) = ChromKMT(i,:);
        end
        if allObj_CD_KMT(i,3) > recordGbest
            recordGbest = allObj_CD_KMT(i,3);
            tourGbest = ChromKMT(i,:);
        end
    end
    
    L_best(N) = recordGbest;
    L_average(N,:) = mean(allObj_CD_KMT, 1);
    fprintf('Iteration %d: recordGbest = %f, tourGbest = %s\n', i, recordGbest, mat2str(tourGbest));
    
end  
