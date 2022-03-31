function accuracy = predict_smoke(J, Y, topk)
% function accuracy = predict_smoke(J, Y, topk)
%
% @param  J         Approximate affinity matrix
% @param  Y         Labels
% @param  topk      k-nearest neighbors prediction (default 1)
%
% @return accuracy  Prediction accuracy
%
% <Reference>
% Yu, Fangchen. "Online Affinity Matrix Calibration", 2022.

if (nargin < 3)
    topk = 1;
end

n = size(J, 1);
J(1:n+1:n*n) = 0;
accuracy = 0;
for i = 1 : n
    v = J(:, i);
    [~, idx] = sort(v, 'descend');
    I = idx(1:topk);
    y_predict = mode(Y(I));
    y_true = Y(i);
    if y_true == y_predict
        accuracy = accuracy + 1;
    end
end
accuracy = accuracy / n;

end