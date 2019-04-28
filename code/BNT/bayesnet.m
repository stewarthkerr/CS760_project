clear all
addpath(genpath("./bnt-master"))

%% Input data
data = csvread("insilico_size10_1_timeseries.csv");
data = data(:, 2:end)';
m = mean(data, 2);
data = data > m;
[~, s] = size(data);
[trainInd, testInd, ~] = dividerand(s, 0.9, 0.1, 0);
trainSet = data(:, trainInd);
testSet = data(:, testInd);

%% Learn simple bayes net
ns = ones(1,10)*2;
[sampled_graphs, accept_ratio, num_edges] = learn_struct_mcmc(trainSet, ns);
[~, idx] = max(accept_ratio);
bnet = mk_bnet(sampled_graphs{idx}, ns);
for i = 1:10
    bnet.CPD{i} = tabular_CPD(bnet, i);
end
learn_params(bnet, trainSet);

%% Test result
testSet = num2cell(testSet);
hide = rand(size(testSet)) > 0.5;
[I,J] = find(hide);
for k=1:length(I)
  testSet{I(k), J(k)} = [];
end
engine = pearl_inf_engine(bnet);
%newengine = enter_evidence(engine, {testSet{:, 2}});
%marginal_nodes(newengine, 0)
%[~, LL] = learn_params_em(engine, testSet, 1);
%find_mpe(engine, {testSet{:, 1}})
%marginal_nodes(engine, nodes)
