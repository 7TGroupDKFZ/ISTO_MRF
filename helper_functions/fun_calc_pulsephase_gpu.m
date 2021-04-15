% function:     fun_calc_pulsephase_gpu
% purpose:      precalculate arrays to give excitation pulse a phase
% inputs:   	number of isochromats, pulse phase
% outputs:      pulse phase arrays, indices (tells when is which phase array needed)

% 01.03.2021 - f.kratzer@dkfz.de

%%

function [phi_out, phi_indices] = fun_calc_pulsephase_gpu(num_of_spins,phi_in)

    % calculate each phase only once
    [phi_unique, ~,phi_indices] = unique(phi_in);

    jmax = length(phi_unique);
    for j=1:jmax
        if phi_unique(j)==0
            phi_out(j).calcphase = 0;
        else
            phi_out(j).calcphase = 1;
            phi = ones([1 num_of_spins])*phi_unique(j); 
            cos_1 = gpuArray(cos(1.*phi));
            cos_2 = gpuArray(cos(2.*phi));
            cos_3 = gpuArray(cos(3.*phi));
            sin_1 = gpuArray(sin(1.*phi));
            sin_2 = gpuArray(sin(2.*phi));
            sin_3 = gpuArray(sin(3.*phi));

            cos_ar = cat(1,ones(size(cos_1)),ones(size(cos_1)),cos_1,cos_1,ones(size(cos_1)),cos_1,cos_1...
                            ,cos_2,cos_2,ones(size(cos_1)),cos_1,cos_1,cos_2,cos_2,cos_3,cos_3);

            sin_ar = cat(1,zeros(size(sin_1)),zeros(size(sin_1)),sin_1,sin_1,zeros(size(sin_1)),sin_1,sin_1...
                            ,sin_2,sin_2,zeros(size(sin_1)),sin_1,sin_1,sin_2,sin_2,sin_3,sin_3);

            phi_out(j).sin_phi = sin_ar;
            phi_out(j).cos_phi = cos_ar;
        end
    end
    clear cos_1 cos_2 cos_3 sin_1 sin_2 sin_3