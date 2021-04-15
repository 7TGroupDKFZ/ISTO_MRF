% function:     fun_relax_J0_matrix
% purpose:      J0 matrix
% inputs:   	J0
% outputs:      J0 matrix

% 01.03.2021 - f.kratzer@dkfz.de

%%
function relax_J0_matrix = fun_relax_J0_matrix(J0)



    relax_J0_matrix = [0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,  0.6*J0        ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    , sqrt(6)*0.2*J0 ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,  0.6*J0        ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           , sqrt(6)*0.2*J0 ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,   J0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,   J0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...    
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,   J0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,   J0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    , sqrt(6)*0.2*J0 ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,  0.4*J0        ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           , sqrt(6)*0.2*J0 ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,  0.4*J0        ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,   J0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,   J0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0;...
                       0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0    ,    0    ,    0    ,    0           ,    0           ,    0    ,    0    ,    0    ,    0];
