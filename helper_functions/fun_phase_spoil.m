% function:     fun_phase_spoil
% purpose:      apply phase spoiling 
% inputs:   	Tlm, precalculated spoiling angles
% outputs:      phase spoiled Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%
function Tlm = fun_phase_spoil(Tlm,cos_spoil,sin_spoil)

    A = Tlm([1 2 4 3 5 7 6 9 8 10 12 11 14 13 16 15],:);
    Tlm = Tlm.*cos_spoil + 1i*A.*sin_spoil;