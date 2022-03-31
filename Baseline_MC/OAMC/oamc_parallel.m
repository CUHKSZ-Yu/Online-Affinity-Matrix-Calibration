function [Jon_cal] = oamc_parallel(Joff, Jon, s, V)
% function [Jon_cal] = oamc_parallel(Joff, Jon)
%
% @param  Joff     Initial offline matrix
% @param  Jon      Initial online matrix
% @param  s        Eigenvalues
% @param  V        Singular matrix
%
% @return Jon_cal  Calibrated online matrix
%
% <Reference>
% Yu, Fangchen. "Online Affinity Matrix Calibration", 2022.

if (nargin < 4)
    [V, S, ~] = svd(Joff);
    s = diag(S);
end

n = size(Jon, 2);
tol = 1e-4;
C = V * diag(sqrt(s));
Cinv = diag(1./s) * C';
U = Cinv * Jon;
U_cal = U;

for i = 1 : n
    u0 = U(:, i);
    U_cal(:, i) = oamc_step(s, u0, tol);
end
Jon_cal = C * U_cal;

end

