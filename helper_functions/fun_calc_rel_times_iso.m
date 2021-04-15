% function:     fun_calc_rel_times_iso
% purpose:      convert spectral density parameters into relaxation parameters with monoexponential T1 (without residual quadrupolar interaction)
% inputs:   	spectral density parameters
% outputs:      relaxation times

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [T1,T2l,T2s] = fun_calc_rel_times_iso(J0, J1, J2);
     
    T2s = 1./(J0+J1);
    T2l = 1./(J1+J2);
    % 1st order Taylor series approximation for monoexponential T1 is used (T1 = 1/(0.8/T1l+0.2/T1s) )
    T1 = 1./(0.8.*(2.*J2) + 0.2.*(2.*J1));