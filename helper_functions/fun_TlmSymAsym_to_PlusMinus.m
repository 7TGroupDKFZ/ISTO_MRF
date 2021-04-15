% function:     fun_TlmSymAsym_to_PlusMinus
% purpose:      transform symmetric and asymmetric tensors (Tlm(s) & Tlm(a)) to Tl+m & Tl-m
% inputs:   	Tlm(s) & Tlm(a)
% outputs:      Tl+m & Tl-m

% 01.03.2021 - f.kratzer@dkfz.de

%%
function [Tlm_pm] = fun_TlmSymAsym_to_PlusMinus(Tlm)

    T00_hist = Tlm(1,:,:);
    T10_hist = Tlm(2,:,:);
    T1m1_hist = (Tlm(3,:,:)+Tlm(4,:,:))*sqrt(2)/2;
    T1p1_hist = (Tlm(3,:,:)-Tlm(4,:,:))*sqrt(2)/2;
    T20_hist = Tlm(5,:,:);
    T2m1_hist = (Tlm(6,:,:)+Tlm(7,:,:))*sqrt(2)/2;
    T2p1_hist = (Tlm(6,:,:)-Tlm(7,:,:))*sqrt(2)/2;
    T2m2_hist = (Tlm(8,:,:)+Tlm(9,:,:))*sqrt(2)/2;
    T2p2_hist = (Tlm(8,:,:)-Tlm(9,:,:))*sqrt(2)/2;
    T30_hist = Tlm(10,:,:);
    T3m1_hist = (Tlm(11,:,:)+Tlm(12,:,:))*sqrt(2)/2;
    T3p1_hist = (Tlm(11,:,:)-Tlm(12,:,:))*sqrt(2)/2;
    T3m2_hist = (Tlm(13,:,:)+Tlm(14,:,:))*sqrt(2)/2;
    T3p2_hist = (Tlm(13,:,:)-Tlm(14,:,:))*sqrt(2)/2;
    T3m3_hist = (Tlm(15,:,:)+Tlm(16,:,:))*sqrt(2)/2;
    T3p3_hist = (Tlm(15,:,:)-Tlm(16,:,:))*sqrt(2)/2;

    Tlm_pm.tensors = cat(1,T00_hist,T10_hist,T1m1_hist,...
                        T1p1_hist,T20_hist,T2m1_hist,T2p1_hist,...
                        T2m2_hist,T2p2_hist,T30_hist,T3m1_hist,...
                        T3p1_hist,T3m2_hist,T3p2_hist,T3m3_hist,T3p3_hist);
    Tlm_pm.titles = {'T00','T10','T1-1',...
                    'T1+1','T20','T2-1','T2+1',...
                    'T2-2','T2+2','T30','T3-1',...
                    'T3+1','T3-2','T3+2','T3-3','T3+3'};