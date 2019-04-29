clear all
addpath(genpath("./GlobalMIT_1.0_Release"))

%% Input data
data = csvread("insilico_size10_1_timeseries.csv");
data = data(:, 2:end);
m = mean(data, 1);
data = double(data > m)+1;
data1 = data(1:21, :);
[a1, b1] = multi_time_series_cat(data1);
data2 = data(22:42, :);
[a2, b2] = multi_time_series_cat(data2);
data3 = data(43:63, :);
[a3, b3] = multi_time_series_cat(data3);
data4 = data(64:84, :);
[a4, b4] = multi_time_series_cat(data4);
data5 = data(85:105, :);
[a5, b5] = multi_time_series_cat(data5);
a = [a1; a2; a3; a4; a5];
b = [b1; b2; b3; b4; b5];
nodeName = num2cell(1:10);

%% Input ground truth
true_net = GroundTruthMat("DREAM4_GoldStandard_InSilico_Size10_1.csv", 10, 10);
s_true=score_MIT(data,true_net);

%% Learn DBN structure
alpha=0.999;
allowSelfLoop=1;
[best_net]=globalMIT_ab(a, b, alpha, allowSelfLoop);
csvwrite("structure.csv", best_net)
