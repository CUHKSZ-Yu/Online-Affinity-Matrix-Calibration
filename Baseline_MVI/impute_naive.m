function [imputedX] = impute_naive(X, method)
% function [imputedX] = impute_naive(X, method)
%
% Replace NaN values by row means/medians/modes.
%
% @param X          Missing dataset
% @param method     Default 'mean' (mean imputation)
% 
% @return imputedX  Imputed matrix with all data samples

if (nargin < 2)
    method = 'mean';
end

idx = isnan(X);
X0 = X;
X0(idx) = 0;
if strcmp(method, 'mean')
    M = nanmean(X, 2);
    imputedX = X0 + M .* idx;
elseif strcmp(method, 'median')
    M = nanmedian(X, 2);
    imputedX = X0 + M .* idx;
elseif strcmp(method, 'mode')
    M = mode(X, 2);
    imputedX = X0 + M .* idx;
end
imputedX(imputedX<=0.5) = 0;
imputedX(imputedX>0.5) = 1;
end