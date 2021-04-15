% function:     fun_rect_pulse_apply
% purpose:      rect pulse
% inputs:   	Tlm, sample parameters, FA, pulse phase, pulse duration
% outputs:      Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [Tlm] = fun_rect_pulse_apply(Tlm,params,FA,phi,pulse_dur)

    % excitation
    params.w1 = FA/pulse_dur;
    % pulse phase
    if phi~=0
        Tlm = fun_rot_Tlm(Tlm,-phi);
    end
    % apply FA pulse with precess function
    [Tlm] = fun_precess(Tlm,params,pulse_dur);
    % pulse phase
    if phi~=0
        Tlm = fun_rot_Tlm(Tlm,phi);
    end