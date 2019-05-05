clear all
addpath(genpath("./GlobalMIT_1.0_Release"))

%% Setting
bin_num = 3;

%% Read data
% Data file
data = readtable(...
    "../../data/Arab.Meristem/arabidopsis.meristem.expression.csv", ...
    'ReadRowNames',true);
data_no0 = data(any(data{:,:},2),:);
% Cluster filte
cluster = readtable("../../output/Clustering.csv", ...
    'ReadRowNames',true);
cluster.Properties.DimensionNames{1} = 'Feature_ID';
cluster.Properties.VariableNames{1} = 'Cluster';
% Combine data ans cluster
data_cluster = join(data_no0, cluster, 'Keys', 'Feature_ID');
final_data = sortrows(varfun(@mean, data_cluster, 'GroupingVariables', 'Cluster'));
cluster_num = size(final_data,1);
% Discritize dataset
final_data_mat = table2array(final_data(:, 3:end));
final_data_dst = zeros(size(final_data_mat));
for i=1:cluster_num
    final_data_dst(i, :) = discretize(final_data_mat(i,:), bin_num);
end
% Seperate dataset
data_1 = final_data_dst(:, 1:2:end);
data_2 = final_data_dst(:, 2:2:end);
% Construct the a & b used for GlobalMIT
[a1, b1] = multi_time_series_cat(data_1');
[a2, b2] = multi_time_series_cat(data_2');
a = [a1; a2];
b = [b1; b2];

%% Learn DBN structure
alpha=0.5;
allowSelfLoop=1;
[best_net]=globalMIT_ab(a, b, alpha, allowSelfLoop);
csvwrite("structure.csv", best_net)

%% Input ground truth
true_net = zeros(cluster_num, cluster_num);
true_index = readtable(...
    "../../data/Arab.Meristem/arabidopsis.meristem.modules.interactions.tsv", ...
    'FileType', 'text', 'Delimiter', '\t', 'HeaderLines', 0);
for cell=table2array(true_index)'
    true_net(int8(cell(1)+1), int8(cell(2)+1)) = cell(3);
end
curr_best_net = best_net & (1-diag(ones(1,cluster_num)));
TP = sum(sum(curr_best_net & true_net));
FP = sum(sum(curr_best_net & ~true_net));
FN = sum(sum(~curr_best_net & true_net));
precission = TP/(TP+FP);
recall = TP/(TP+FN);
acc = 1-(FP+FN)/(cluster_num*(cluster_num-1));

%% Precision recall curve
% PR = zeros(100,2);
% for i=1:100
%     alpha = 1-exp(-i*0.1);
%     [curr_best_net] = globalMIT_ab(a, b, alpha, allowSelfLoop);
%     curr_best_net = curr_best_net & (1-diag(ones(1,cluster_num)));
%     TP = sum(sum(curr_best_net & true_net));
%     FP = sum(sum(curr_best_net & ~true_net));
%     FN = sum(sum(~curr_best_net & true_net));
%     precission = TP/(TP+FP);
%     recall = TP/(TP+FN);
%     acc = 1-(FP+FN)/(cluster_num*(cluster_num-1));
%     PR(i, 1) = alpha;
%     PR(i, 2) = precission;
%     PR(i, 3) = recall;
%     PR(i, 4) = acc;
% end
% plot(PR(:,3), PR(:,2))