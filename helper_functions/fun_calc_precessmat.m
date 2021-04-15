% function:     fun_calc_precessmat
% purpose:      precalculate matrices needed for precession
% inputs:   	params with spectral density parameters
% outputs:      precalculated matrices S and Md, needed for precession

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [S,Md] = fun_calc_precessmat(params)

    % w1 corresponds to excitation
    if~isfield(params,'w1')
        params.w1 = 0;
    end
    [S,Md] = eig(fun_excite_matrix(params.w1) ...
            + fun_offres_matrix(params.woff) ...
            + fun_quadrupolar_matrix(params.wq) ...
            - fun_relax_J0_matrix(params.J0) ...
            - fun_relax_J1J2_matrix(params.J1, params.J2));

