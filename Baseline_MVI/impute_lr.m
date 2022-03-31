function [imputedX, time] = impute_lr(Xoff, Xon)
% function [imputedX, time] = impute_lr(Xoff, Xon)
%
% Impute a data matrix. Each column is a sample. Each NaN value in a vector 
% is replaced by the statistical value calculated by Linear Regression from
% observed values to missing values. (see reference)
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% 
% @return imputedX  Imputed matrix with all data samples
% @return time      Running time
%
% <Reference>
% Seber, George AF, and Alan J. Lee. Linear regression analysis. John Wiley & Sons, 2012.

n_off = size(Xoff, 2);
n_on = size(Xon, 2);
imputedX = [Xoff, Xon];
if mod(n_on, 100) == 0
    time = zeros(1, n_on/100);
end

tic;
for i = 1 : n_on
    y = Xon(:,i);
    idx = isnan(y);
    y_obs = y(~idx);
    X_mis = Xoff(idx,:);
    X_obs = Xoff(~idx,:);
    beta = regress(y_obs, [ones(sum(~idx),1), X_obs]);
    y_mis = sum(beta' .* [ones(sum(idx),1), X_mis], 2);
    y(idx) = y_mis;
    imputedX(:, n_off+i) = y;
    if mod(i, 100) == 0
        time(1, i/100) = toc;
    end
end
imputedX(imputedX<=0.5) = 0;
imputedX(imputedX>0.5) = 1;
end