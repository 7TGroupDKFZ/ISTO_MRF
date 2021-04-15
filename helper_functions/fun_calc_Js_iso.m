% function:     fun_calc_Js_iso
% purpose:      convert T1,T2l and T2s in spectral denisity parameters (without residual quadrupolar interaction)
% inputs:   	relaxation times [s]
% outputs:      spectral density parameters

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [J0 J1 J2] = fun_calc_Js_iso(T1,T2l,T2s)

    % 1st order Taylor series approximation for monoexponential T1 is used (T1 = 1/(0.8/T1l+0.2/T1s) )
     J0 = 5./(6.*T1) - 4./(3.*T2l) + 1./T2s;
     J1 = 4./(3.*T2l) - 5./(6.*T1);
     J2 = 5./(6.*T1) - 1./(3.*T2l);
