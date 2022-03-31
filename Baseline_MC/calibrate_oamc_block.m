function [Joamc_block] = calibrate_oamc_block(Jmiss, noff, non, model, maxiter)
% function [Joamc_block] = calibrate_oamc_block(Jmiss, noff, non, model, maxiter)
%
% @param Jmiss          Initial affinity matrix
% @param noff           Number of offline samples
% @param non            Number of online samples
% @param model          Type of batch calibration method ('DMC' or 'CMC')
% @param maxiter        Maximum iterations
%
% @return Joamc_block   Calibrated affinity matrix
%
% <Reference>
% Yu, Fangchen. "Online Affinity Matrix Calibration", 2022.

if (nargin < 5)
    maxiter = 50;
elseif (nargin < 4)
    model = 'cmc';
end

Jon = Jmiss(noff+1:end, noff+1:end);
if strcmp(model, 'dmc')
    Jon_cal = calibrate_dmc(Jon, maxiter);
elseif strcmp(model, 'cmc')
    Jon_cal = calibrate_cmc(Jon, maxiter);
end

Joff = Jmiss(1:noff, 1:noff);
Jpar_ini = Jmiss(1:noff, noff+1:end);
Jpar_cal = oamc_parallel(Joff, Jpar_ini);

Joamc_block = [Jmiss(1:noff, 1:noff), Jpar_cal; Jpar_cal', Jon_cal];

end


