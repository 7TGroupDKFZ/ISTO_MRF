% function:     fun_rot_Tlm
% purpose:      rotate array
% inputs:   	Tlm, rotation angle phi [rad]
% outputs:      Tlm

% 01.03.2021 - f.kratzer@dkfz.de

%%
function Tlm = fun_rot_Tlm(Tlm,phi)

    Tlm_array_temp = Tlm;

    sin_phi1 = sin(1.*phi);
    sin_phi2 = sin(2.*phi);
    sin_phi3 = sin(3.*phi);

    cos_phi1 = cos(1.*phi);
    cos_phi2 = cos(2.*phi);
    cos_phi3 = cos(3.*phi);

% 	Tlm_array_temp( 1,:) = Tlm( 1,:); % T00 unaffected, only here for completeness
% 	Tlm_array_temp( 2,:) = Tlm( 2,:); % T10 unaffected, only here for completeness
    Tlm_array_temp( 3,:) = Tlm( 3,:).*cos_phi1 + 1i.*Tlm( 4,:).*sin_phi1;
    Tlm_array_temp( 4,:) = Tlm( 4,:).*cos_phi1 + 1i.*Tlm( 3,:).*sin_phi1;
%   Tlm_array_temp( 5,:) = Tlm( 5,:); % T20 unaffected, only here for completeness
    Tlm_array_temp( 6,:) = Tlm( 6,:).*cos_phi1 + 1i.*Tlm( 7,:).*sin_phi1;
    Tlm_array_temp( 7,:) = Tlm( 7,:).*cos_phi1 + 1i.*Tlm( 6,:).*sin_phi1;
    Tlm_array_temp( 8,:) = Tlm( 8,:).*cos_phi2 + 1i.*Tlm( 9,:).*sin_phi2;
    Tlm_array_temp( 9,:) = Tlm( 9,:).*cos_phi2 + 1i.*Tlm( 8,:).*sin_phi2;
%   Tlm_array_temp(10,:) = Tlm(10,:); % T30 unaffected, only here for completeness
    Tlm_array_temp(11,:) = Tlm(11,:).*cos_phi1 + 1i.*Tlm(12,:).*sin_phi1;
    Tlm_array_temp(12,:) = Tlm(12,:).*cos_phi1 + 1i.*Tlm(11,:).*sin_phi1;
    Tlm_array_temp(13,:) = Tlm(13,:).*cos_phi2 + 1i.*Tlm(14,:).*sin_phi2;
    Tlm_array_temp(14,:) = Tlm(14,:).*cos_phi2 + 1i.*Tlm(13,:).*sin_phi2;
    Tlm_array_temp(15,:) = Tlm(15,:).*cos_phi3 + 1i.*Tlm(16,:).*sin_phi3;
    Tlm_array_temp(16,:) = Tlm(16,:).*cos_phi3 + 1i.*Tlm(15,:).*sin_phi3;

    Tlm = Tlm_array_temp;