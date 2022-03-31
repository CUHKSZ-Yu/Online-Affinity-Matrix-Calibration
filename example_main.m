%% 

clear all; clc; warning off;
addpath(genpath('Baseline_MVI'));
addpath(genpath('Baseline_MC'));

load('data_demo.mat');
dataset = 'Data_demo'; 

%% Parameter Setting
p_list = [0.3, 0.5, 0.7]; % values of missing ratio p
niter = 5;                % number of iterations
noff = 1000;              % number of offline samples
non = 100;                % number of online samples
seed = 2022;              % random seed for reproducing results

fprintf('\nMLJ special issue 2022 submission "Online Affinity Matrix Calibration"');
fprintf('\nDemo: missing data processing in Section 5.3\n');

for k = 1 : length(p_list)
    p = p_list(k);
    fprintf(['\n',dataset]); fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f', noff, non, p, niter);
    rng(seed);
    
    %% Online Scenario
    for i = 1 : niter
        fprintf('\nIters = %1.0f: ', i);
        
        %% Data Construction
        n = size(X, 2);
        rp = randperm(n);
        Xoff = double(full(X(:, rp(1:noff))));
        Xon = double(full(X(:, rp(noff+1:noff+non))));
        Xtrue = [Xoff, Xon];
        Jtrue = jaccard_approx(Xtrue); 

        Xon(rand(size(Xon)) < p) = NaN;  
        Xmiss = [Xoff, Xon];
        Jmiss = jaccard_approx([Xoff, Xon], 'miss'); 
        Fnorm = norm(Jmiss-Jtrue, 'fro')^2;
        Abs = sum(sum(abs(Jmiss-Jtrue)));
        
        %% Online Process
        % ====================== ZERO Imputation ==========================
        fprintf('ZERO, ');
        tic; Xzero = impute_zero(Xmiss); time(i,1) = toc;
        Jzero = jaccard_approx(Xzero);
        rmse(i,1) = norm(Jzero-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,1) = sum(sum(abs(Jzero-Jtrue))) / Abs;
        
        % ====================== MEAN Imputation ==========================
        fprintf('MEAN, ');
        [Xmean, Tmean] = impute_mean(Xoff, Xon); time(i,2) = Tmean(end);
        Jmean = jaccard_approx(Xmean);
        rmse(i,2) = norm(Jmean-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,2) = sum(sum(abs(Jmean-Jtrue))) / Abs;

        % ====================== kNN Imputation ===========================
        fprintf('kNN, '); 
        [Xknn, Tknn] = impute_knn(Xoff, Xon, Jmiss); time(i,3) = Tknn(end);
        Jknn = jaccard_approx(Xknn);
        rmse(i,3) = norm(Jknn-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,3) = sum(sum(abs(Jknn-Jtrue))) / Abs;
        
        % ====================== LR Imputation ============================
        fprintf('LR, '); 
        [Xlr, Tlr] = impute_lr(Xoff, Xon); time(i,4) = Tlr(end);
        Jlr = jaccard_approx(Xlr);
        rmse(i,4) = norm(Jlr-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,4) = sum(sum(abs(Jlr-Jtrue))) / Abs;
        
        % ====================== GROUSE Imputation ========================
        fprintf('GROUSE, '); 
        [Xgr, Tgr] = impute_grouse(Xmiss); time(i,5) = Tgr(end);
        Jgr = jaccard_approx(Xgr);
        rmse(i,5) = norm(Jgr-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,5) = sum(sum(abs(Jgr-Jtrue))) / Abs;
        
        % ====================== KFMC Imputation ==========================
        fprintf('KFMC, '); 
        tic; Xkfmc = impute_kfmc(Xmiss, 'on', 'rbf', Xtrue); time(i,6) = toc;
        Jkfmc = jaccard_approx(Xkfmc);
        rmse(i,6) = norm(Jkfmc-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,6) = sum(sum(abs(Jkfmc-Jtrue))) / Abs;
        
        % ====================== DMC Calibration ==========================
        fprintf('DMC, '); 
        tic; Jdmc = calibrate_dmc(Jmiss); time(i,7) = toc;
        rmse(i,7) = norm(Jdmc-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,7) = sum(sum(abs(Jdmc-Jtrue))) / Abs;

        % ====================== CMC Calibration ==========================
        fprintf('CMC, '); 
        tic; Jcmc = calibrate_cmc(Jmiss); time(i,8) = toc;
        rmse(i,8) = norm(Jcmc-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,8) = sum(sum(abs(Jcmc-Jtrue))) / Abs;
        
        % ====================== OAMC Calibration =========================
        fprintf('OAMC, '); 
        [Joamc, Toamc] = calibrate_oamc(Jmiss, noff, non); time(i,9) = Toamc(end);
        rmse(i,9) = norm(Joamc-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,9) = sum(sum(abs(Joamc-Jtrue))) / Abs;
        
        % ====================== OAMC-DMC Calibration =====================
        fprintf('OAMC_DMC, '); 
        tic; Joamcd = calibrate_oamc_block(Jmiss, noff, non, 'dmc'); time(i,10) = toc;
        rmse(i,10) = norm(Joamcd-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,10) = sum(sum(abs(Joamcd-Jtrue))) / Abs;
        
        % ====================== OAMC-CMC Calibration =====================
        fprintf('OAMC_CMC, '); 
        tic; Joamcc = calibrate_oamc_block(Jmiss, noff, non, 'cmc'); time(i,11) = toc;
        rmse(i,11) = norm(Joamcc-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,11) = sum(sum(abs(Joamcc-Jtrue))) / Abs;
        
        fprintf('Finish.');
    end

    %%
    fprintf(['\n',dataset]); fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f\n', noff, non, p, niter);

    stat = [mean(rmse); mean(rmae); mean(time)];
    Stat = roundn(stat, -4);
    Table = table(Stat(:,1),Stat(:,2),Stat(:,3),Stat(:,4),Stat(:,5),Stat(:,6),Stat(:,7),Stat(:,8),Stat(:,9),Stat(:,10),Stat(:,11),...
         'VariableNames',{'ZERO','MEAN','kNN','LR','GROUSE','KFMC','DMC','CMC','OAMC','OAMC_DMC','OAMC_CMC'},...
         'RowNames',{'RMSE';'RMAE';'Time'});
    disp(Table)
end



