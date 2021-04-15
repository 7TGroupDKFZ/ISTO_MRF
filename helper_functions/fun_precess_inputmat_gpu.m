% function:     fun_precess_inputmat_gpu
% purpose:      spin evolution with or without presence of pulse; precalculated S and Md
% inputs:   	Tlm, evolution duration, pulse duration, S, Md
% outputs:      Tlm 

% 01.03.2021 - f.kratzer@dkfz.de

%%

function Tlm = fun_precess_inputmat_gpu(Tlm_in,duration,S,Md)
     
    A = pagefun(@mtimes,S,pagefun(@mtimes, pagefun(@times,eye(16),exp(duration*Md)), pagefun(@inv,S)));
    Tlm = pagefun(@mtimes,A,Tlm_in);
