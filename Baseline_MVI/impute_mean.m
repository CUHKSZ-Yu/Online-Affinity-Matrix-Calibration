function [imputedX, time] = impute_mean(Xoff, Xon)
% function [imputedX, time]] = impute_mean(Xoff, Xon)
%
% Replace NaN values by row means.
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% 
% @return imputedX  Imputed matrix with all data samples
% @return time      Running time

n_off = size(Xoff, 2);
n_on = size(Xon, 2);
if mod(n_on, 100) == 0
    time = zeros(1, n_on/100);
end

tic;
imputedX = [Xoff, Xon];
rowmean = mean(Xoff, 2);
for i = 1 : n_on
    x = Xon(:, i);
    idx = isnan(x);
    imputedX(idx, n_off+i) = rowmean(idx);
    if mod(i, 100) == 0
        time(1, i/100) = toc;
    end
end
imputedX(imputedX<=0.5) = 0;
imputedX(imputedX>0.5) = 1;
end