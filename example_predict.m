%% 

clear all; clc; warning off; 
addpath(genpath('Baseline_MVI'));
addpath(genpath('Baseline_MC'));

load('data_predict.mat');
dataset = 'Data_predict';

%% Density Cutoff
n = size(X,2);
nan_idx = isnan(X);
nan_num = sum(nan_idx,2);
density = 1 - nan_num/n;
X0 = X;
idx = find(density>=0.3);
X3 = X(idx,:);
idx = find(density>=0.5);
X5 = X(idx,:);
idx = find(density>=0.7);
X7 = X(idx,:);
X_list{1} = X0; X_list{2} = X3; X_list{3} = X5; X_list{4} = X7;

%% Parameter Setting
density_list = [0, 0.3, 0.5, 0.7];
noff = 2000; non = n - noff;
fprintf('\nMLJ special issue 2022 submission "Online Affinity Matrix Calibration"');
fprintf('\nDemo: prediction task in Section 6.3\n');

%% Prediction Task
for j = 1:length(density_list)
    fprintf('\nIters = %1.0f: ', j);
    Xmiss = X_list{j};
    density = density_list(j);
    Jmiss = jaccard_approx(Xmiss, 'miss');
    accuracy(j,1) = predict_smoke(Jmiss, Y);
    
    % ====================== ZERO Imputation ==========================
    fprintf('ZERO, ');
    Xzero = impute_zero(Xmiss);
    Jzero = jaccard_approx(Xzero);
    accuracy(j,2) = predict_smoke(Jzero, Y);
    
    % ====================== MEAN Imputation ==========================
    fprintf('MEAN, ');
    Xmean = impute_naive(Xmiss, 'mean');
    Jmean = jaccard_approx(Xmean);
    accuracy(j,3) = predict_smoke(Jmean, Y);

    % ====================== kNN Imputation ===========================
    fprintf('kNN, '); 
    Xknn = knnimpute(Xmiss, 10); Xknn(Xknn<=0.5) = 0; Xknn(Xknn>0.5) = 1;
    Jknn = jaccard_approx(Xknn);
    accuracy(j,4) = predict_smoke(Jknn, Y);

    % ====================== GROUSE Imputation ========================
    fprintf('GROUSE, '); 
    [Xgr, ~] = impute_grouse(Xmiss); 
    Jgr = jaccard_approx(Xgr);
    accuracy(j,5) = predict_smoke(Jgr, Y);

    % ====================== KFMC Imputation ==========================
    fprintf('KFMC, '); 
    Xkfmc = impute_kfmc(Xmiss, 'off', 'poly');
    Jkfmc = jaccard_approx(Xkfmc);
    accuracy(j,6) = predict_smoke(Jkfmc, Y);
    
    % ====================== DMC Calibration ==========================
    fprintf('DMC, ');
    Jdmc = calibrate_dmc(Jmiss); 
    accuracy(j,7) = predict_smoke(Jdmc, Y);
    Jmiss(1:noff, 1:noff) = Jdmc(1:noff, 1:noff);

    % ====================== OAMC Calibration =========================
    fprintf('OAMC, '); 
    [Joamc, ~] = calibrate_oamc(Jmiss, noff, non);
    accuracy(j,8) = predict_smoke(Joamc, Y);

    % ====================== OAMC-DMC Calibration =====================
    fprintf('OAMC-DMC, '); 
    Joamcd = calibrate_oamc_block(Jmiss, noff, non, 'dmc', 5);
    accuracy(j,9) = predict_smoke(Joamcd, Y);

    fprintf('finished.');
end

%%
fprintf('\nPrediction accuracy of smoking behavior on NHANES dataset:\n');
Stat = roundn(accuracy, -4);
Table = table(Stat(:,1),Stat(:,2),Stat(:,3),Stat(:,4),Stat(:,5),Stat(:,6),Stat(:,7),Stat(:,8),Stat(:,9),...
        'VariableNames',{'J^0','ZERO','MEAN','kNN','GROUSE','KFMC','DMC','OAMC','OAMC-DMC'},...
        'RowNames',{'cutoff=0.0';'cutoff=0.3';'cutoff=0.5';'cutoff=0.7'});
disp(Table)


