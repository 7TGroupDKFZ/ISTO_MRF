% function:     fun_rot_Tlm_gpu
% purpose:      rotate array
% inputs:   	Tlm, rotation angle phi [rad]
% outputs:      Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%

function Tlm = fun_rot_Tlm_gpu(cos_phi,sin_phi,Tlm)

	Tlm = Tlm.*cos_phi + 1i* Tlm([1 2 4 3 5 7 6 9 8 10 12 11 14 13 16 15],:,:).*sin_phi;
        
