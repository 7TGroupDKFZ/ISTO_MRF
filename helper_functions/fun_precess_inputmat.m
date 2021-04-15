% function:     fun_precess_inputmat
% purpose:      spin evolution with or without presence of pulse; precalculated S and Md 
% inputs:   	Tlm, evolution duration, pulse duration, S, Md
% outputs:      Tlm 


% 01.03.2021 - f.kratzer@dkfz.de

%%
function Tlm = fun_precess_inputmat(Tlm_array,duration,S,Md)

    Tlm = (S* (eye(16).*exp(duration*Md)) *inv(S)) * Tlm_array;





