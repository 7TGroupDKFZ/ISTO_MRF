% function:     fun_Lorentzian_distribution
% purpose:      apply Lorentzian distribution in B0-direction to introduce T2*
% inputs:   	T2-signal, relaxation parameters, B0, minT2l_fact: limits Lorentzian (T2l>T2l*>T2l * minT2l_fact)
% outputs:      T2*-signal, lookup tables with: 1:biexponential T1, 2:monoexponential T1, 3: in J-basis

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [signal_Lorentz, LUT_T_biT1_Lorentz, LUT_T_monT1_Lorentz,LUT_J_Lorentz] = fun_Lorentzian_distribution(signal,LUT_Parameter,minT2l_fact);
    
    

    % entries that are supposed to appear in the dictionary
    B0_range_calc = -50:2:50; %min(B0_range_calc)>=! min(B0)-100Hz so that the Lorentz distribution fits

    % B0 interpolation step size relative to B0 input 
    B0_stepsize = 0.5; %1Hz input -> 0.5Hz sampling

    % B0 range evaluated for Lorentzian distribution
    x = (-100:B0_stepsize:100); % Lorentzian evaluated from -100Hz+B0_center to +100+B0_center Hz
    %% calculate T1,T2l,T2s and B0

    % convert J0, J1, J2 into relaxation times
    [T1l, T1s,T2l,T2s] = fun_calc_rel_times_iso_biT1(LUT_Parameter(:,1), LUT_Parameter(:,2), LUT_Parameter(:,3));
    B0 = LUT_Parameter(:,5)/2/pi;
    B0_u = unique(B0);

    signal_temp = permute(reshape(signal,size(signal,1),length(B0_u),[]),[1 3 2]);
    clear signal_cat

    %% interpolate B0

    steps = 1/B0_stepsize;
    signal_temp_interp = complex(single(zeros([size(signal_temp,1) size(signal_temp,2) (size(signal_temp,3)-1)*steps])));
    % linear interpolation along B0-axis
    for j = 1:size(signal_temp,2)
        for i = 1:(size(signal_temp,3)-1)
           signal_temp_interp(:,j,((i-1)*steps + (1:steps) )) =  permute(1/(steps)*((steps:-1:1)'*transpose(signal_temp(:,j,i))+(0:(steps-1))'*transpose(signal_temp(:,j,i+1))),[2 3 1]);
           B0_u2(((i-1)*steps + (1:steps) )) = 1/(steps)*((steps:-1:1)*B0_u(i) + (0:(steps-1))*B0_u(i+1));
        end
    end
    signal_temp_interp(:,:,end+1) = single(signal_temp(:,:,end));
    B0_u2(end+1) = B0_u(end);
    clear signal_temp;
    
    % needed for convolution
    B0_start = find(round(B0_u2,2) == round(B0_range_calc(1),2));
    % setup B0 distribution and reshape signal
    signal_temp_interp = permute(signal_temp_interp,[1 3 2]);


    %% initialize the new signal evolutions

    % count the number of entries for initialization
    signal_size = 0;
    for i_cnt = 1:size(signal_temp_interp,3)
        T2l_set = T2l(i_cnt);
        kmax = (T2l_set*1e3-round(T2l_set*1e3*minT2l_fact));
        for k_cnt = 1:kmax
            signal_size = signal_size+1;
        end
    end
    signal_size = signal_size*length(B0_range_calc);

    % final signal (uncompressed)
    signal_Lorentz = complex(single(zeros([size(signal_temp_interp,1) signal_size])));
    % lookup table with biexponential T1: T1l, T1s, T2l*, T2s*,B0
    LUT_T_biT1_Lorentz = single(zeros([signal_size 5]));
    % lookup table with monoexponential T1:  T1, T2l*, T2s*,B0
    LUT_T_monT1_Lorentz = single(zeros([signal_size 4]));
    % lookup table in J-basis: J0,J1,J2,B0, width of Lorentzian
    LUT_J_Lorentz = single(zeros([signal_size 5]));

    %% run combination
    l_cnt = 1;
    i_cnt_max = size(signal_temp_interp,3);
    fprintf('Lorentz distribution progress: %i%% ',0)
    for i_cnt = 1:i_cnt_max
        T1l_set = T1l(1+(i_cnt-1)*length(B0_u));
        T1s_set = T1s(1+(i_cnt-1)*length(B0_u));
        T2l_set = T2l(1+(i_cnt-1)*length(B0_u));
        T2s_set = T2s(1+(i_cnt-1)*length(B0_u));
        J0_set = LUT_Parameter(1+(i_cnt-1)*length(B0_u),1);
        J1_set = LUT_Parameter(1+(i_cnt-1)*length(B0_u),2);
        J2_set = LUT_Parameter(1+(i_cnt-1)*length(B0_u),3);


        kmax = (T2l_set*1e3-round(T2l_set*1e3*minT2l_fact));
        for k_cnt = 1:kmax
            % calculate lorentzian distribution
            T2l_Star = round(T2l_set*minT2l_fact+(k_cnt-1)*1e-3,3);
            g = (1/T2l_Star - 1/T2l_set) / (2*pi);
            f = 1/(pi*g) * (g^2./(x.^2+g^2));
            f = f/sum(f);
            T2s_Star = round((g*2*pi + 1/(T2s_set))^-1,4);
            temp = conv2(f,signal_temp_interp(:,:,i_cnt));
            % save signal for specified B0_range_calc
            signal_Lorentz(:,l_cnt:l_cnt+length(B0_range_calc)-1) = temp(:,floor(length(f)/2)+B0_start:(B0_range_calc(2)-B0_range_calc(1))/B0_stepsize:floor(length(f)/2)+B0_start+(length(B0_range_calc)-1)*(B0_range_calc(2)-B0_range_calc(1))/B0_stepsize);
            % lookup tables
            LUT_J_Lorentz(l_cnt:l_cnt+length(B0_range_calc)-1,:) = cat(2,repmat([J0_set J1_set J2_set],length(B0_range_calc),1),B0_range_calc.',repmat(g,length(B0_range_calc),1));
            LUT_T_biT1_Lorentz(l_cnt:l_cnt+length(B0_range_calc)-1,:) = cat(2,repmat([T1l_set*1e3 T1s_set*1e3 T2l_Star*1e3 T2s_Star*1e3],length(B0_range_calc),1),B0_range_calc.');
            LUT_T_monT1_Lorentz(l_cnt:l_cnt+length(B0_range_calc)-1,:) = cat(2,repmat([1./(0.8./T1l_set+0.2./T1s_set)*1e3 T2l_Star*1e3 T2s_Star*1e3],length(B0_range_calc),1),B0_range_calc.');
            l_cnt = l_cnt+length(B0_range_calc);
        end
        % show distribution calculation status
        if sum(floor(i_cnt_max/10)*(1:9)-i_cnt == 0) == 1
            fprintf('%i%% ',10*i_cnt/floor(i_cnt_max/10))
        end
    end
    fprintf('%i%% \n',100)

