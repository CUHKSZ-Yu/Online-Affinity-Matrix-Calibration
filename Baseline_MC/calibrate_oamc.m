function [J, time] = calibrate_oamc(Jmiss, noff, non)
% function [J, time] = calibrate_oamc(Jmiss, noff, non)
%
% @param Jmiss  Initial affinity matrix
% @param noff   Number of offline samples
% @param non    Number of online samples
%
% @return J     Calibrated affinity matrix
%
% <Reference>
% Yu, Fangchen. "Online Affinity Matrix Calibration", 2022.

if mod(non, 100) == 0
    time = zeros(1, non/100);
end

Jobs = Jmiss(1:noff, 1:noff);
Jimp = Jobs;

tic;
for i = 1 : non
    v = Jmiss(1:noff+i-1, noff+i);
    vonl = onlcal(Jimp, v);
    Jimp = [Jimp, vonl; vonl', 1];
    if mod(i, 100) == 0
        time(1, i/100) = toc;
    end
end
J = Jimp;

end


%% One-step OAMC
function [vonl] = onlcal(J0, v0)
% function [vonl] = onlcal(J0, v0)
%
% @param  S0    Pairwise similarity matrix
% @param  v0    Intial similarity vector
% @return vonl  Calibrated similarity vector

n = size(J0, 1);
tol = 1e-4;

[U, S, V] = svd(J0);
s = diag(S);
C = U * diag(sqrt(s));
Cinv = diag(1./s) * C';

y0 = Cinv * v0;
if norm(y0) <= 1
    ycal = y0;
else
    lambda_min = max(sqrt(s'.^2 * y0.^2) - max(s), 0);
    lambda_max = sqrt(s'.^2 * y0.^2) - min(s);
    lambda = lambda_min;
    ycal = (s ./ (s+lambda)) .* y0;
    ylen = norm(ycal);
    while (ylen > 1) || (ylen < 1-tol)
        lambda = 0.5*(lambda_min + lambda_max);
        ycal = (s ./ (s+lambda)) .* y0;
        ylen = norm(ycal);
        % use bisection search to find optimal lambda
        if ylen > 1
            lambda_min = lambda;
        elseif ylen < 1-tol
            lambda_max = lambda;
        end 
    end
end
vonl = C * ycal;
end
