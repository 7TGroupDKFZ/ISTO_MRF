% function:     fun_calc_excitation_gpu
% purpose:      precalculate the excitation matrices, used in fun_rect_pulse_apply_gpu
% inputs:   	sample parameters, pulse duration, FA array
% outputs:      excitation matrix, FA indices (tells when is which excitation matrix needed)

% 01.03.2021 - f.kratzer@dkfz.de

%%

function [ex_mat, FA_indices] = fun_calc_excitation_gpu(params, pulse_dur, FA)

    % calculate each pulse only once
    [FA_unique, ~,FA_indices] = unique(FA);
    w1 = FA_unique /pulse_dur;

    % init matrices
    ex_mat = zeros([16 16 length(params.woff) length(FA_unique)]);
    jmax = length(FA_unique);
    kmax = length(params.woff);
    % calculate excitation matrix to be used in fun_rect_pulse_apply_gpu
    for k = 1:kmax
        for j = 1:jmax
            [S,Md] = eig(fun_excite_matrix(w1(j)) ...
                    + fun_offres_matrix(params.woff(k)) ...
                    + fun_quadrupolar_matrix(params.wq(k)) ... 
                    - fun_relax_J0_matrix(params.J0(k)) ...
                    - fun_relax_J1J2_matrix(params.J1(k), params.J2(k)));
            ex_mat(:,:,k,j) = ((S)* (diag(ones([16,1])).*exp(pulse_dur*Md)) *inv(S));
        end
    end
