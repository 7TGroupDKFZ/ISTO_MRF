% function:     fun_precalc_spoiling_gpu
% purpose:      precalculate spoiling angles, can be used for spoiling with 'fun_rot_Tlm_gpu'
% inputs:   	number of isochromats
% outputs:      spoiling angles

% 01.03.2021 - f.kratzer@dkfz.de

%%

function [cos_spoil, sin_spoil] = fun_precalc_spoiling_gpu(num_of_spins)

    phi = -pi:2*pi/num_of_spins:(pi-2*pi/num_of_spins); 
    cos_spoil_1 = gpuArray(cos(1.*phi));
    cos_spoil_2 = gpuArray(cos(2.*phi));
    cos_spoil_3 = gpuArray(cos(3.*phi));
    sin_spoil_1 = gpuArray(sin(1.*phi));
    sin_spoil_2 = gpuArray(sin(2.*phi));
    sin_spoil_3 = gpuArray(sin(3.*phi));
    
    cos_spoil = cat(1,ones(size(cos_spoil_1)),ones(size(cos_spoil_1)),...
                cos_spoil_1,cos_spoil_1,ones(size(cos_spoil_1)),cos_spoil_1,...
                cos_spoil_1,cos_spoil_2,cos_spoil_2,ones(size(cos_spoil_1)),...
                cos_spoil_1,cos_spoil_1,cos_spoil_2,cos_spoil_2,cos_spoil_3,cos_spoil_3);

    sin_spoil = cat(1,zeros(size(sin_spoil_1)),zeros(size(sin_spoil_1)),...
                sin_spoil_1,sin_spoil_1,zeros(size(sin_spoil_1)),sin_spoil_1,...
                sin_spoil_1,sin_spoil_2,sin_spoil_2,zeros(size(sin_spoil_1)),...
                sin_spoil_1,sin_spoil_1,sin_spoil_2,sin_spoil_2,sin_spoil_3,sin_spoil_3);
    clear cos_spoil_1 cos_spoil_2 cos_spoil_3 sin_spoil_1 sin_spoil_2 sin_spoil_3