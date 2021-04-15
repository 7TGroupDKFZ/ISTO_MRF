% function:     fun_rect_pulse
% purpose:      rect pulse
% inputs:   	Tlm, sample parameters, FA [rad], pulse phase [rad], pulse duration [s], pulse samples, pulse acquisition (is sampling of Tlm during pulse wanted or not)
% outputs:      Tlm, acquired Tlm history, time steps

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [Tlm, Tlm_history, t_pulse] = fun_rect_pulse(Tlm,params, FA,phi,pulse_dur,pulse_samples,pulse_acq)

    % if Tlm is to be recorded during pulse
    if pulse_acq == 1
        Tlm_history = complex(zeros([size(Tlm,1),1,pulse_samples]));
        t_pulse = zeros([1 pulse_samples]);
        for i_pulse = 1:pulse_samples
            Tlm = fun_rect_pulse_apply(Tlm,params, FA/pulse_samples,phi,pulse_dur/pulse_samples );
            Tlm_history(:,1,i_pulse) = sum(Tlm,2);
            t_pulse(i_pulse) = i_pulse*pulse_dur/pulse_samples;

        end
    % if Tlm is not to be recorded, Tlm at the end of the pulse can be calculated directly
    else
        Tlm_history = [];
        t_pulse = [];
        pulse_samples = 1;
        Tlm = fun_rect_pulse_apply(Tlm,params, FA/pulse_samples,phi,pulse_dur/pulse_samples);
    end


    