% function:     MRF_simulation
% purpose:      simulate MRF fingerprints (unaccelerated: no GPU or parallel computing in here) 
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
%% set parameters

% example data (calculate (J0,J1,J2) = (180 40 20)Hz and (220 30 15)Hz for off-resonances: -150 to 150 Hz)
J0 = cat(1,ones([301 1])*180, ones([301 1])*220);         	%[Hz]
J1 = cat(1,ones([301 1])*40, ones([301 1])*30);          	%[Hz]
J2 = cat(1,ones([301 1])*20, ones([301 1])*15);             %[Hz]
wq = cat(1,ones([301 1])*0, ones([301 1])*0);               %[Hz]
woff = cat(1,(-150:1:150).',(-150:1:150).')*2*pi;           %off-resonance [Hz]*2*pi

% cat parameters
LUT_Parameter = cat(2,J0,J1,J2,wq, woff);

num_of_spins   = 100;                                       %number of isochromats
%% set FA and TE pattern
pulse_duration = 1e-3;          % pulse duration [s]
pulse_phase = zeros([1e3 1]);   % pulse phase [rad]
TE = 0.55;                      % minimal echso time [s]
B1 = 1.0;                       % relative B1 amplitude

% TE pattern
% TE_Var corresponds to TE - pulse_duration/2
load('helper_functions/TE_Var_Na23')  
TE_scal = 20;
TE_Var = (ceil((TE_Var*TE_scal)*100)/100+TE- pulse_duration/2*1e3)*1e-3; % time between pulse and acquisition: TE_scal*TE_Var + TE- d_pulse/2:

% TE_b and double echo pattern
load('helper_functions/dE');
% TE_Var_b corresponds to TE_b - pulse_duration/2
TE_Var_b = TE_Var_b*1e-3;

% FA pattern
load('helper_functions/FA_Pattern');
FA_scal = 90;
FA_Pattern = FA_Pattern*FA_scal/100/180*pi*B1; %scale FA pattern

TR_fil = 13.79*1e-3; % time between sampling of k_space center and next pulse [s]

%% initialize  

% precalculate phase spoiling
[cos_spoil, sin_spoil] = fun_precalc_spoiling(num_of_spins);

% signal
num_entries = size(LUT_Parameter,1);
signal = complex(zeros([signal_length num_entries]));

%% run simulation
fprintf('simulation progress: %i%% ',0)
for i_num = 1:num_entries
    % init
    Tlm = zeros([16 num_of_spins]);
    Tlm(1,:) = 1;
    Tlm(2,:) = 1; 
    params.J0 = LUT_Parameter(i_num,1);
    params.J1 = LUT_Parameter(i_num,2);
    params.J2 = LUT_Parameter(i_num,3);
    params.wq = LUT_Parameter(i_num,4);
    params.woff = LUT_Parameter(i_num,5);
    
    [S,Md] = fun_calc_precessmat(params);    % precalculate precession matrices
    
    cycle_ind = 1;
    % loop over all cycles
    for FA_ind =1:length(FA_Pattern)
        % excitation
        Tlm = fun_rect_pulse((Tlm),params, FA_Pattern(FA_ind),pulse_phase(FA_ind),pulse_duration,1,0);
        % double echo?
        if double_echo(FA_ind)==1
            % TEb
            Tlm = fun_precess_inputmat(Tlm,TE_Var_b(FA_ind),S,Md);
            % ADC
            [signal_p, Tlm] = fun_acquire_point(Tlm);
            signal(cycle_ind,i_num) = signal_p;
            cycle_ind = cycle_ind+1;
            % TE
            Tlm = fun_precess_inputmat(Tlm,TE_Var(FA_ind)-TE_Var_b(FA_ind),S,Md);
            % ADC
            [signal_p, Tlm] = fun_acquire_point(Tlm);
            signal(cycle_ind,i_num) = signal_p;
            cycle_ind = cycle_ind+1;
        else
            % TE
            Tlm = fun_precess_inputmat(Tlm,TE_Var(FA_ind),S,Md);
            % ADC
            [signal_p, Tlm] = fun_acquire_point(Tlm);
            signal(cycle_ind,i_num) = signal_p;
            cycle_ind = cycle_ind+1;
        end
        % TR_fil
        Tlm = fun_precess_inputmat(Tlm,TR_fil,S,Md);
        % spoiling
        Tlm = fun_phase_spoil(Tlm,cos_spoil,sin_spoil);
        
    end
    
	% show simulation status
	if sum(floor(num_entries/10)*(1:9)-i_num == 0) == 1
        fprintf('%i%% ',10*i_num/floor(num_entries/10))
	end
        
end
fprintf('%i%% \n',100)

%clear temp variables
clearvars -except LUT_Parameter signal
%% Lorentz distribution on signal to introduce T2*

minT2l_fact = 0.6;  % T2l>T2l*>T2l * minT2l_fact

% in this function (fun_Lorentzian_distribution) a B0-range of -150:1:150Hz is assumed and has to be adapted if this is changed
[signal_Lorentz, LUT_T_biT1_Lorentz, LUT_T_monT1_Lorentz,LUT_J_Lorentz] = fun_Lorentzian_distribution(signal,LUT_Parameter,minT2l_fact);

