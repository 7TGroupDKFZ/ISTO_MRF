% function:     fun_calc_rel_times_iso_biT1
% purpose:      convert spectral density parameters into relaxation parameters with biexponential T1 (without residual quadrupolar interaction)
% inputs:   	spectral density parameters
% outputs:      relaxation times

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [T1l, T1s,T2l,T2s] = fun_calc_rel_times_iso_biT1(J0, J1, J2);
 T2s = 1./(J0+J1);
 T2l = 1./(J1+J2);
 T1l = 1./(2.*J2); 
 T1s = 1./(2.*J1);