% function:     fun_acquire_point_gpu
% purpose:      acquire current signal
% inputs:   	Tlm
% outputs:      signal

% 01.03.2021 - f.kratzer@dkfz.de

%%

function [signal] = fun_acquire_point_gpu(Tlm3,Tlm4)

    signal = gather(sum(Tlm4-Tlm3,2));

                    

                        
