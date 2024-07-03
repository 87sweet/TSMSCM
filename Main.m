clear;clc;format compact;tic;hold on

% 生成数据范围
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

% 生成前3列和后9列的行数
front_rows = 30;
back_rows = 30;

aa=3;
bb=9;

% 生成符合正态分布的随机数据
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

% 合并前3列和后9列的数据
Co = [Co_front, Co_back];
D1 = [D1_front, D1_back];  %F1
D2 = [D2_front, D2_back];  %F2
D3 = [D3_front, D3_back];  %F3
T = [T_front, T_back];  %T
C = [C_cost_front, C_cost_back];  %C
Re = [Re_front, Re_back];  %Re
Av = [Av_front, Av_back];  %Av
Sta = [Sta_front, Sta_back];  %Sta

% 进行最大最小归一化处理
%min_max_normalize = @(x) (x - min(x(:))) / (max(x(:)) - min(x(:)));
max_min_normalize = @(x) (max(x(:))-x) / (max(x(:)) - min(x(:)));

T_normalized = max_min_normalize(T);
C_normalized = max_min_normalize(C);

% 保存数据
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

% 保存为.mat文件
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

%% 参数设置
NIND=20;        %种群大小
NINDO=20;        %种群大小 
MAXGEN=100;     %迭代次数
Pc=0.9;         %交叉概率
Pm=0.09;        %变异概率
%N=size(C,2);

%% 初始化种群 

%% 优化
N=1;
ChromKMT_xnew1=InitPopKMT(NIND,C_KMT);
xnew1=ChromKMT_xnew1;
allObj_CD_KMT_xnew1=allObject_CD_KMT(ChromKMT_xnew1,C_KMT,D1_KMT,D2_KMT,D3_KMT); 
[value,index]=max(allObj_CD_KMT_xnew1(:,3)); 
tourPbest=ChromKMT_xnew1; 
tourGbest=ChromKMT_xnew1(index,:); 
recordPbest=-inf*ones(NIND,1);               %个体最优记录
recordGbest=allObj_CD_KMT_xnew1(index,3);                %群体最优记录
L_best=zeros(MAXGEN,1);          %记录每次迭代过程中全局最优个体的总距离  %L_best(N)=allObj_CD_KMT(index1);
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
    recordPbestO=-inf*ones(NINDO,1);               %个体最优记录
    recordGbestO=allObj_EF_OMT_xnew2(index111,5);                %群体最优记录
    L_bestO=zeros(MAXGEN,1);          %记录每次迭代过程中全局最优个体的总距离  %L_best(N)=allObj_CD_KMT(index1);
    L_averageO=zeros(MAXGEN,size(allObj_EF_OMT_xnew2, 2));
    
    ChromOMT=InitPopOMT(NINDO,E1_OMT);
    allObj_EF_OMT=allObject_EF_OMT(ChromOMT,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back); 
    

    for M=1:MAXGEN
        for i=1:NINDO
            allObj_EF_KMT=allObject_EF_KMT(ChromKMT,E1_KMT,E2_KMT,E3_KMT,E4_KMT,F_KMT,T_front,C_cost_front);
            
            a = 250; % 示例值，请根据实际情况设置/时间约束
            b = 16000; % 示例值，请根据实际情况设置/成本约束

            % 获取 allObj_EF_OMT 和 allObj_EF_KMT 第三列的值
            OMT_col_time = allObj_EF_OMT(:, 1);
            KMT_col_time = allObj_EF_KMT(:, 1);
            total_time=OMT_col_time+KMT_col_time;
            OMT_col_cost = allObj_EF_OMT(:, 2);
            KMT_col_cost = allObj_EF_KMT(:, 2);
            total_cost=OMT_col_cost+ KMT_col_cost;
            % 检查每行的值是否满足约束条件
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
                       
             % 1111与个体最优进行交叉
            L=size(E1_OMT,2);
            c1=randi([1 L],1,1); %产生交叉位
            c2=randi([1 L],1,1); %产生交叉位
            while c1==c2
                c2=randi([1 L],1,1);
            end
            
                        % 定义和复制初始数组元素
            a0 = xnew2(i,:);  % 复制xnew2的第i行到a0，用作备份和参考
            b0 = tourPbestO(i,:);  % 复制tourPbestO的第i行到b0

            % 确定交叉点的边界
            chb1 = min(c1, c2);  % 确定交叉点的下界
            chb2 = max(c1, c2);  % 确定交叉点的上界

            % 交叉循环
            for j = chb1:chb2
                % 在xnew2的第i行第j列放入b0的第j个元素
                xnew2(i,j) = b0(j);  
            end

            % 查找并修正可能的重复元素
            for j = chb1:chb2
                x = find(xnew2(i,:) == xnew2(i,j));  % 查找xnew2第i行中与当前交叉位置元素相等的所有位置
                i1 = x(x ~= j);  % 从找到的位置中排除当前交叉点的位置
                if ~isempty(i1)  % 如果存在重复的元素
                    % 用原始的a0中的元素替换掉重复元素，确保每个元素是唯一的
                    xnew2(i,i1) = a0(j);  
                end
            end
            
            dist=allObject_EF_OMT(xnew2,E1_OMT,E2_OMT,E3_OMT,E4_OMT,F_OMT,T_back,C_cost_back); %allObjectKMT 求目标的函数
            if  allObj_EF_OMT(i,5)<dist(i,5) %allObjKMT为目标函数allObjectKMT求到的目标值
            ChromOMT(i,:)=xnew2(i,:);
            end 
            
            % 交叉操作
            c1 = randi([1 L], 1, 1); % 产生交叉位
            c2 = randi([1 L], 1, 1); % 产生交叉位

            if c1 ~= c2
                a0 = xnew2(i, :);  % 保存原始数据行的状态
                b0 = tourGbestO;   % 获取群体最优解

                chb1 = min(c1, c2);  % 确定交叉点的下界
                chb2 = max(c1, c2);  % 确定交叉点的上界

                % 进行交叉操作
                for j = chb1:chb2
                    xnew2(i, j) = b0(j);  % 将b0的第j个元素复制到xnew2的第i行第j列
                end

                % 查找并修正可能的重复元素
                for j = chb1:chb2
                    x = find(xnew2(i, :) == xnew2(i, j));  % 查找xnew2第i行中与当前交叉位置元素相等的所有位置
                    i1 = x(x ~= j);  % 从找到的位置中排除当前交叉点的位置
                    if ~isempty(i1)  % 如果存在重复的元素
                        xnew2(i, i1) = a0(j);  % 用原始的a0中的元素替换掉重复元素
                    end
                end
            end

            % 计算新解的目标函数值并比较
            dist2 = allObject_EF_OMT(xnew2, E1_OMT, E2_OMT, E3_OMT, E4_OMT, F_OMT, T_back, C_cost_back);
            if allObj_EF_OMT(i) < dist2
                ChromOMT(i, :) = xnew2(i, :);
            end
            
        
        %% 变异操作
            c1=randi([1 L],1,1); %产生交叉位
            c2=randi([1 L],1,1); %产生交叉位
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
            % 获取 allObj_EF_OMT 和 allObj_EF_KMT 第三列的值
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

    
%% 将下层模型最优解输入上层模型，求解上层模型的最优解
    for i=1:NIND
        allObj_CD_KMT=allObject_CD_KMT(ChromKMT,C_KMT,D1_KMT,D2_KMT,D3_KMT);  
        allObj_CD_OMT=allObject_CD_KMT(tourGbestO,C_OMT,D1_OMT,D2_OMT,D3_OMT);       
        c = 10.6; % 示例值，请根据实际情况设置/时间约束

        % 获取 allObj_EF_OMT 和 allObj_EF_KMT 第三列的值
        OMT_col_Co = allObj_CD_OMT(:,1);
        KMT_col_Co = allObj_CD_KMT(:,1);
        total_Co=allObj_CD_KMT(:, 1)+allObj_CD_OMT(:,1);
        % 检查每行的值是否满足条件
        for i = 1:length(OMT_col_Co)
            if total_Co(i)< c
                allObj_CD_OMT(i, 3) = 0;
            end
        end
        
        % 1111与个体最优进行交叉
        L=size(C_KMT,2);
        c1=randi([1 L],1,1); %产生交叉位
        c2=randi([1 L],1,1); %产生交叉位
        while c1==c2
            c2=randi([1 L],1,1);
        end
            
         % 定义和复制初始数组元素
         a0 = xnew1(i,:);  % 复制xnew2的第i行到a0，用作备份和参考
         b0 = tourPbest(i,:);  % 复制tourPbestO的第i行到b0

         % 确定交叉点的边界
         chb1 = min(c1, c2);  % 确定交叉点的下界
         chb2 = max(c1, c2);  % 确定交叉点的上界

         % 交叉循环
         for j = chb1:chb2
             % 在xnew2的第i行第j列放入b0的第j个元素
              xnew1(i,j) = b0(j);  
         end

            % 查找并修正可能的重复元素
         for j = chb1:chb2
              x = find(xnew1(i,:) == xnew1(i,j));  % 查找xnew2第i行中与当前交叉位置元素相等的所有位置
              i1 = x(x ~= j);  % 从找到的位置中排除当前交叉点的位置
              if ~isempty(i1)  % 如果存在重复的元素
                    % 用原始的a0中的元素替换掉重复元素，确保每个元素是唯一的
                  xnew1(i,i1) = a0(j);  
              end
          end
            
          dist=allObject_CD_KMT(xnew1,C_KMT,D1_KMT,D2_KMT,D3_KMT); %allObjectKMT 求目标的函数
          if  allObj_CD_KMT(i,3)<dist(i,3) %allObjKMT为目标函数allObjectKMT求到的目标值
          ChromKMT(i,:)=xnew1(i,:);
          end 

       % 2222与群体最优进行交叉
        c1=randi([1 L],1,1); %产生交叉位
        c2=randi([1 L],1,1); %产生交叉位
        if c1~=c2
            a0=xnew1(i,:);b0=tourGbest; 
            chb1=min(c1,c2);     
            chb2=max(c1,c2);
            for j=chb1:chb2
                a1=xnew1(i,:); b1=tourGbest;
                xnew1(i,j)=b0(j);  
                x=find(xnew1(i,:)==xnew1(i,j));  %查找第i行中与a(j)值相等的元素的位置，并将结果保存在向量x中。
                i1=x(x~=j); %从向量x中排除位置为j的元素，将结果保存在向量i1中。
                if ~isempty(i1)  %判断向量i1是否为空。
                xnew1(i,i1)=a1(j); %如果向量i1非空，将向量a中位置为i1的元素更新为向量a1中的第i个元素。
                end
            end
        end
        dist2=0;     
        dist2=dist2+allObject_CD_KMT(xnew1,C_KMT,D1_KMT,D2_KMT,D3_KMT);
        if  allObj_CD_KMT<dist2
            ChromKMT(i,:)=xnew1(i,:);
        end
        
        %% 变异操作
        c1=randi([1 L],1,1); %产生交叉位
        c2=randi([1 L],1,1); %产生交叉位
        while c1==c2
%             c1=randi([1 L],1,1); %产生交叉位
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
        c = 11; % 示例值，请根据实际情况设置/时间约束

        % 获取 allObj_EF_OMT 和 allObj_EF_KMT 第三列的值
        OMT_col_Co = allObj_CD_OMT(:,1);
        KMT_col_Co = allObj_CD_KMT(:,1);
        total_Co=allObj_CD_KMT(:, 1)+allObj_CD_OMT(:,1);
        % 检查每行的值是否满足条件
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
