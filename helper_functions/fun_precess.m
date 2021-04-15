% function:     fun_precess
% purpose:      spin evolution with or without presence of pulse
% inputs:   	Tlm, sample parameters, pulse duration
% outputs:      Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [Tlm] = fun_precess(Tlm,params,pulse_dur)

    % excitation
    if~isfield(params,'w1')
        params.w1 = 0;
    end
    [S,Md] = eig(fun_excite_matrix(params.w1) ...
            + fun_offres_matrix(params.woff) ...
            + fun_quadrupolar_matrix(params.wq) ...
            - fun_relax_J0_matrix(params.J0) ...
            - fun_relax_J1J2_matrix(params.J1, params.J2));
    Tlm = (S* (eye(16).*exp(pulse_dur*Md)) *inv(S))* Tlm;





