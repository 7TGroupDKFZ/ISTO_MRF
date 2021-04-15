% function:     fun_init_Tlm_gpu
% purpose:      initiate Tlm
% inputs:   	number of spins, number of parameter combinations
% outputs:      Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%

function Tlm = fun_init_Tlm_gpu(fov_spins, num_comb);

    % initialize Tlm 
    Tlm_o = gpuArray(ones([1 fov_spins]));
    Tlm_z = gpuArray(zeros([1 fov_spins]));
    Tlm = single(cat(1,Tlm_o,Tlm_o,Tlm_z,Tlm_z,Tlm_z,Tlm_z,Tlm_z,Tlm_z...
                 ,Tlm_z,Tlm_z,Tlm_z,Tlm_z,Tlm_z,Tlm_z,Tlm_z,Tlm_z));
    Tlm = repmat(Tlm,1,1,num_comb);