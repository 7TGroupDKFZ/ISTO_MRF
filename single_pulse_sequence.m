% function:     single_pulse_sequence
% purpose:      easy framework example to simulate an FID
%
% The theory is based on: Hancu I, Van der Maarel J, Boada F. A Model for the Dynamics of Spins 3/2 in Biological Media:Signal Loss during Radiofrequency Excitation in Triple-Quantum-Filtered Sodium MRI. Journal of Magnetic Resonance 2000; 147(2):179-191.
%
% ISTO simulations use symmetric and assymetric Tlm: Tlm(s) and Tlm(a) as described by Hancu et al.
%
% The simulations are based on the work by Arthur Magill: https://github.com/arthurmagill/isto_sim
%
% Fabian Kratzer - German Cancer Research Center (DKFZ)
% 01.03.2021 - f.kratzer@dkfz.de

%%
addpath('helper_functions')

%% set parameters:

% J-parameters can also be calculated via: [J0 J1 J2] = fun_calc_Js_iso(T1,T2l,T2s) - T1,T2l,T2s in [s]

% sampel parameters
params.J0 =180;             %[Hz]
params.J1 = 20;             %[Hz]
params.J2 = 20;             %[Hz]
params.wq = 0;              %[Hz]
params.woff = 0*2*pi;       %off-resonance [Hz]*2*pi

num_of_spins   = 1;         %number of isochromats

pulse_FA = 90 * pi/180;     % FA [rad]
pulse_phase = 0;            % pulse phase
pulse_duration = 5e-3;      % pulse duration [s]
pulse_samples = 1000;       % pulse samples
pulse_acq = 1;              % acquire Tlm during pulse 1/0

acq_duration = 50e-3;       % acquisition duration after pulse [s]
acq_dwelltime = 1.0e-4;     % acquisition dwelltime [s]
%% initialize tensors
Tlm = zeros([16 num_of_spins]);
Tlm(1,:) = 1;
Tlm(2,:) = 1;
%% FID sequence:
% rect pulse
[Tlm, Tlm_history_pulse,t_pulse] = fun_rect_pulse(Tlm,params, pulse_FA,pulse_phase,pulse_duration,pulse_samples,pulse_acq);

% precession 
[Tlm_history_prec, Tlm,t_prec] = fun_precess_inputmat_acquire(Tlm,params,acq_dwelltime,acq_duration/acq_dwelltime);
t_prec = t_prec+pulse_duration;

% concatenate all timings and signal evolutions:
t = cat(2,t_pulse,t_prec);
Tlm_history = cat(3,Tlm_history_pulse,Tlm_history_prec);



%% plot results

% transform Tlm(a) and Tlm(s) to Tl+m and Tl-m (just for plotting)
Tlm_hist_pm = fun_TlmSymAsym_to_PlusMinus(Tlm_history);


figure
for i_plot = 1:16
    subplot(4,4,i_plot);
    plot(t(:)*1e3,squeeze(round(real(sum(Tlm_hist_pm.tensors(i_plot,:,:),2)),15)));
    title(Tlm_hist_pm.titles{i_plot})
    hold on
    plot(t(:)*1e3,squeeze(round(imag(sum(Tlm_hist_pm.tensors(i_plot,:,:),2)),15)));
    xlabel('t (ms)')
    hold off; box on;
    maxlim = max(abs(ylim));
    ylim([-maxlim maxlim]);
    %line where the pulse ends
    line([t_pulse(end)*1e3 t_pulse(end)*1e3],[-maxlim maxlim],'color','k')
    legend('real', 'imag')
end
sgtitle('spherical tensor operators')


figure
plot(t_prec*1e3,-squeeze(round(real(sum(Tlm_hist_pm.tensors(4,:,end+1-length(t_prec):end),2)),15))*2/sqrt(2));
hold on
plot(t_prec*1e3,-squeeze(round(imag(sum(Tlm_hist_pm.tensors(4,:,end+1-length(t_prec):end),2)),15))*2/sqrt(2));
legend('real', 'imag')
hold off; box on;
xlabel('t (ms)')
title('Signal')

