function J = jaccard_approx(X, choice)
% function [J] = jaccard_approx(X, choice))
%
% Approximate a Jaccard index matrix for samples with NaN values.
%
% @param X       d*n, each column is a sample
% @param choice  Default "true" (calculate the true Jaccard index matrix)
% @return J      n*n, Jaccard index matrix

if (nargin < 2)
    choice = 'true';
end

if strcmp(choice, 'miss')
    [d, n] = size(X);
    O = ~isnan(X); % {0,1} matrix: 1 means known feature
    OP = (X==1);   % {0,1} matrix: 1 means observed feature
    P = zeros(n);
    Q = zeros(n);
    for i = 1 : n
        A = repmat(OP(:,i), 1, n);
        P(i,:) = full(sum(A & OP));
        Q(i,:) = full(sum(A & O));
    end
    J = P ./ (Q + Q' - P);
    J(1:n+1:n*n) = 1;
    J(isnan(J)) = 0;
elseif strcmp(choice, 'true')
    [d, n] = size(X);
    P = full(X' * X);
    nf = full(sum(X));
    Q = repmat(nf, n, 1);
    S = Q + Q' - P;
    J = P ./ S;
    J(1:n+1:n*n) = 1;
    J(isnan(J)) = 0;
elseif strcmp(choice, 'basic') % common features
    [d, n] = size(X);
    Idx = isnan(X);
    Xzero = X; Xzero(Idx) = 0;
    P = Xzero' * Xzero;
    for i = 1 : n
        idx = Idx(:,i);
        if sum(idx) == 0
            Q(i,:) = sum(Xzero);
        else
            idx = repmat(idx, 1, n);
            XI = Xzero; XI(idx) = 0;
            Q(i,:) = sum(XI);
        end
    end
    J = P ./ (Q + Q' - P);
    J(1:n+1:n*n) = 1;
    J(isnan(J)) = 0;
end

end