% function:     fun_rect_pulse_apply_gpu
% purpose:      rect pulse
% inputs:   	Tlm, excitation matrices, pulse phase struct
% outputs:      Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%

function [Tlm] = fun_rect_pulse_apply_gpu(Tlm, ex_mat, phi)

    % pulse phase
    if phi.calcphase~=0
        num_of_spins = size(Tlm,2);
        Tlm = fun_rot_Tlm_gpu(phi.cos_phi,-phi.sin_phi,Tlm);
    end
        Tlm =  pagefun(@mtimes, ex_mat, Tlm);
    % pulse phase
    if phi.calcphase~=0
        Tlm = fun_rot_Tlm_gpu(phi.cos_phi,phi.sin_phi,Tlm);
    end