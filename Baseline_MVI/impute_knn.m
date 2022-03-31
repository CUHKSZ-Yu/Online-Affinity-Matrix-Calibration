function [imputedX, time] = impute_knn(Xoff, Xon, Smiss, k)
% function [imputedX, time] = impute_knn(Xoff, Xon, Smiss, k)
%
% Impute a data matrix. Each column is a sample. Each NaN value is replaced
% by the mean of the sample's k-nearest neighbors with known values. If all 
% k-nearest samples' corresponding features are NaN, then replaced by zero.
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% @param k          Default 1 (the top-1 nearest neighbor)
% 
% @return imputedX  Imputed matrix with all data samples
% @return time      Running time

if (nargin < 4)
    k = 1;
end

n_off = size(Xoff, 2);
n_on = size(Xon, 2);
imputedX = [Xoff, Xon];
if mod(n_on, 100) == 0
    time = zeros(1, n_on/100);
end

tic;
nan_idx = isnan(Xon);
for i = 1 : n_on
    s = Smiss(1:n_off, n_off+i);
    [~, idx] = sort(s, 'descend');
    imputedX(nan_idx(:,i), n_off+i) = mean(Xoff(nan_idx(:,i), idx(1:k)), 2);
    if mod(i, 100) == 0
        time(1, i/100) = toc;
    end
end
imputedX(imputedX<=0.5) = 0;
imputedX(imputedX>0.5) = 1;
end