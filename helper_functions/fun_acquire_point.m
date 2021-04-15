% function:     fun_acquire_point
% purpose:      acquire current signal
% inputs:   	Tlm
% outputs:      signal, Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%

function [signal, Tlm] = fun_acquire_point(Tlm)

    signal = sum(Tlm(4,:)-Tlm(3,:),2);
     

                        
