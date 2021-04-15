% function:     fun_calc_precessmat_gpu
% purpose:      precalculate the precession matrices
% inputs:   	sample parameters
% outputs:      S, Md

% 01.03.2021 - f.kratzer@dkfz.de

%%

function [S,Md] = fun_calc_precessmat_gpu(woff,wq,J0,J1,J2)
    
    % initialize
    S = (zeros([16 16 length(woff)]));
    Md = (zeros([16 16 length(woff)]));
    % calculate S and Md for all entries
    for k = 1:length(woff)
        [S(:,:,k),Md(:,:,k)] = eig(fun_offres_matrix(woff(k)) ...
                                + fun_quadrupolar_matrix(wq(k)) ...
                                - fun_relax_J0_matrix(J0(k)) ...
                                - fun_relax_J1J2_matrix(J1(k), J2(k)));
    end