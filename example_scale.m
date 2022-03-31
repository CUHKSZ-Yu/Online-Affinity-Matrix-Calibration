%% 

clear all; clc;
addpath(genpath('Baseline_MC'));

load('data_demo.mat');
dataset = 'Data_demo'; 

%% Parameter Setting
p_list = [0.3, 0.5, 0.7]; % values of missing ratio p
niter = 5;                % number of iterations
noff = 5000;              % number of offline samples
non = 5000;               % number of online samples
seed = 2022;              % random seed for reproducing results

fprintf('\nMLJ special issue 2022 submission "Online Affinity Matrix Calibration"');
fprintf('\nDemo: scalability analysis in Section 5.4\n');

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
        Xoff = X(:, rp(1:noff));
        Xon = X(:, rp(noff+1 : noff+non));
        Xtrue = [Xoff, Xon];
        Jtrue = jaccard_approx(Xtrue);

        Xon(rand(size(Xon)) < p) = NaN;  
        Xmiss = [Xoff, Xon];
        Jmiss = jaccard_approx([Xoff, Xon], 'miss'); 
        Fnorm = norm(Jmiss-Jtrue, 'fro')^2;
        Abs = sum(sum(abs(Jmiss-Jtrue)));
        
        %% SOAMC Calibration
        % ====================== SOAMC-DMC Calibration ====================
        fprintf('SOAMC_DMC, '); 
        tic; Jsoamc = calibrate_soamc(Jmiss, noff, non, 'dmc'); time(i,1) = toc;
        rmse(i,1) = norm(Jsoamc-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,1) = sum(sum(abs(Jsoamc-Jtrue))) / Abs;
        
        % ====================== SOAMC-CMC Calibration ====================
        fprintf('SOAMC_CMC, '); 
        tic; Jsoamc = calibrate_soamc(Jmiss, noff, non, 'cmc'); time(i,2) = toc;
        rmse(i,2) = norm(Jsoamc-Jtrue, 'fro')^2 / Fnorm;
        rmae(i,2) = sum(sum(abs(Jsoamc-Jtrue))) / Abs;
        
        fprintf('Finish.');
    end

    %%
    fprintf(['\n',dataset]); fprintf(': noff=%1.0f, non=%1.0f, p=%1.2f, niter=%1.0f\n', noff, non, p, niter);

    stat = [mean(rmse); mean(rmae); mean(time)];
    Stat = roundn(stat, -4);
    Table = table(Stat(:,1),Stat(:,2),...
         'VariableNames',{'SOAMC_DMC','SOAMC_CMC'},...
         'RowNames',{'RMSE';'RMAE';'Time'});
    disp(Table)
end



