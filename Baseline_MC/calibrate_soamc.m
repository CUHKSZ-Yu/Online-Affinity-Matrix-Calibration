function [Jsoamc] = calibrate_soamc(Jmiss, noff, non, model, maxiter)
% function [Jsoamc] = calibrate_soamc(Jmiss, noff, non, model, maxiter)
%
% @param Jmiss     Initial affinity matrix
% @param noff      Number of offline samples
% @param non       Number of online samples
% @param model     Type of direct calibration method ('DMC' or 'CMC')
% @param maxiter   Maximum iterations (default 10)
%
% @return Jsoamc   Calibrated affinity matrix
%
% <Reference>
% Yu, Fangchen. "Online Affinity Matrix Calibration", 2022.

if (nargin < 5)
    maxiter = 10;
elseif (nargin < 4)
    model = 'dmc';
end
koff = 1000; % size of offline submatrix (default 1000)
kon = 1000;  % size of online submatrix (default 1000)

npart_off = noff / koff;
Joff_par = Jmiss(1:noff, noff+1:end);
for i = 1 : npart_off
    scale = [(i-1)*koff+1 : i*koff];
    Joff = Jmiss(scale, scale);
    Jpar_ini = Jmiss(scale, noff+1:end);
    Jpar_cal = oamc_parallel(Joff, Jpar_ini);
    Joff_par(scale, :) = Jpar_cal;
end

Jon = Jmiss(noff+1:end, noff+1:end);
Jon_cal = calibrate_on_para(Jon, kon, model, maxiter);

Jsoamc = [Jmiss(1:noff, 1:noff), Joff_par; Joff_par', Jon_cal];

end