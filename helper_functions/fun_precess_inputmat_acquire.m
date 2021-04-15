% function:     fun_precess_inputmat_acquire
% purpose:      Tlm evolution is saved under free precession
% inputs:   	Tlm, sampling rate, number of samples
% outputs:      acquired Tlm history, final Tlm, time steps corresponding to Tlm history 

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [Tlm_history, Tlm,t_prec] = fun_precess_inputmat_acquire(Tlm,params,t_dwell,acq_pts)

    [S,Md] = eig( fun_offres_matrix(params.woff) ...
            + fun_quadrupolar_matrix(params.wq) ...
            - fun_relax_J0_matrix(params.J0) ...
            - fun_relax_J1J2_matrix(params.J1, params.J2));
    ev = (S* (diag(ones([16,1])).*exp(t_dwell*Md)) *inv(S));
    signal = complex(zeros([16, 1, acq_pts]));

    % acquire Tlm for all time steps
    for p=1:acq_pts
        Tlm = ev*Tlm;
        Tlm_history(:,:,p) = sum(Tlm,2);
    end

    % time steps 
    t_prec = t_dwell:t_dwell:t_dwell*acq_pts;



                    

                        
