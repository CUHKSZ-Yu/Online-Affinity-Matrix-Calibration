function [J] = calibrate_on_para(Jon, kon, model, maxiter)
% function [J] = calibrate_on_para(Jon, kon, model, maxiter)
%
% @param Jon       Initial online matrix
% @param kon       size of online submatrix (default 1000)
% @param model     Type of direct calibration method ('DMC' or 'CMC')
% @param maxiter   Maximum iterations
%
% @return J        Calibrated online matrix
%
% <Reference>
% Yu, Fangchen. "Online Affinity Matrix Calibration", 2022.

if (nargin < 4)
    maxiter = 50;
elseif (nargin < 3)
    model = 'dmc';
elseif (nargin < 2)
    kon = 1000;
end

non = size(Jon, 1);
npart_on = non / kon;
J = zeros(non, non);
for i = 1 : npart_on
    scale = [(i-1)*kon+1 : i*kon];
    Jon_ini = Jon(scale, scale);
    Jpar_ini = Jon(scale, i*kon+1:non);
    if strcmp(model, 'dmc')
        Jon_cal = calibrate_dmc(Jon_ini, maxiter);
    elseif strcmp(model, 'cmc')
        Jon_cal = calibrate_cmc(Jon_ini, maxiter);
    end
    Jpar_cal = oamc_parallel(Jon_cal, Jpar_ini);
    J(scale, (i-1)*kon+1:end) = [Jon_ini, Jpar_cal];
    J(i*kon+1:end, scale) = Jpar_cal';
end

end


