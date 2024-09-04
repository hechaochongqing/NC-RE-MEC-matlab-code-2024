%clear

load('NO281day2_LFP.mat', 'CSC5_003', 'CSC5_003_ind', 'CSC5_003_ts', 'CSC5_003_ts_step')   %加载LFP的数据
begin_time = CSC5_003_ts(1);   %数据起始时间点
start_time = [4757.07293        %NREM start time
4845.87421
4913.37629
5021.83837
5343.62013
6306.67165
6637.30077
6794.58717
6850.94813
7053.78205

];       %NREM end time
end_time = [4825.88573                
4900.92445
4987.75965
5068.04125
5391.13373
6345.33789
6739.20925
6822.11229
6883.38845
7099.98493

];      
freq = 1/CSC5_003_ts_step;


data_start = zeros(0);
new_freq = 100;    %下采样到new_freq Hz
for t = 1 : length(start_time)
    
    data_temp = CSC5_003((start_time(t)-begin_time)*freq : (end_time(t)-begin_time)*freq);   %1段NREM LFP数据    
    EEG = pop_importdata('dataformat','array','nbchan',1,'data',data_temp','srate', freq);
    EEG.setname='test';
    data_resample = pop_resample(EEG, new_freq);
    data_resample = data_resample.data;
    data = cat(2, data_start, data_resample);
    data_start = data;
    
end

[wt, f] = cwt(data, "amor", new_freq);
wt = abs(wt);
col_theta = find(f>4 & f <12);
col_delta = find(f<4);
mean1 = mean(wt(col_theta, :), 'all');
temp1 = wt(col_theta, :);
std1 = std(temp1(:));

window_size1 = 0.1 * new_freq;   %滑动窗大小（s）
per_all = zeros(length(start_time), 5);  %第一列为epoch时长占比，第二列为epoch数量，第三列为epoch平均时长，第四列为epoch平均强度比例，第五列为非epoch期强度占比
for j = 1 : length(start_time)
    
    bas = 0;
    if j > 1
        for ib = 1 : j-1
            bas = bas + (end_time(ib) - start_time(ib))*new_freq;
        end
    end
    
    window_num1 = floor((end_time(j)-start_time(j))/0.1);
    per_temp = zeros(length(window_num1),2);
    for i = 1 : window_num1
        
        theta_temp = wt(col_theta, bas+window_size1*(i-1)+1:bas+window_size1*i);
        delta_temp = wt(col_delta, bas+window_size1*(i-1)+1:bas+window_size1*i);
        mean2 = mean(theta_temp, 'all');
        if mean2 > mean1 + std1
            per_temp(i,1) = 1;
        end
        
        per_temp(i,2) = sum(theta_temp, 'all')/(sum(theta_temp, 'all') + sum(delta_temp, 'all'));
        
    end
    per_all(j, 1) = sum(per_temp(:,1))/length(per_temp);
    
    col_0 = find(per_temp(:,1)==0);
    per_all(j, 5) = mean(per_temp(col_0, 2));
    
    col_1 = find( per_temp(:,1)==1 );   %寻找连续窗，定义为1个epoch
    c1 = 1;
    arrest = cell(0,0);
    while(c1 < numel(col_1))
        c2 = 0;
        while (c1+c2+1 <= numel(col_1) && col_1(c1)+c2+1 == col_1(c1+c2+1))
           c2 = c2 + 1;
        end
        
        if (c2 >=1)
            arrest = [arrest;(col_1(c1:1:c1+c2))];
        end
        c1 = c1 + c2 +1;
    end
    
    if exist('arrest') && numel(arrest) > 0
        epoch_per = zeros(numel(arrest),2);
        for ii = 1 : numel(arrest)
            epoch_per(ii,1) = numel(arrest{ii,1});
            epoch_per(ii,2) = mean(per_temp(arrest{ii,1},2));
        end
    end
    
    per_all(j, 2) = length(epoch_per);
    per_all(j, 3) = (sum(epoch_per(:,1))*0.1)/per_all(j, 2);
    per_all(j, 4) = sum(epoch_per(:,2))/per_all(j, 2);
    
end

Per_all_L=per_all(:,3:5);


    