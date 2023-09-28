%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%       "Evolutionary dynamics of glucose-deprived cancer cells:      %%%
%%%     insights from experimentally-informed mathematical modelling"   %%%
%%%                                                                     %%%
%%% L. Almeida, J. Denis, N. Ferrand, T. Lorenzi, A. Prunet, M. Sabbah, %%%
%%%                         C. Villa (*)                                %%%
%%%                                                                     %%%
%%%                 Submitted for publication (2023)                    %%%
%%%                  (Preprint available on ArXiv)                      %%%
%%%                                                                     %%%
%%%                                                                     %%%
%%% (*) Email: chiara.villa.1[at]sorbonne-universite.fr                 %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                       %
%%% Almeidaetal2023evolutionary.m: calibration of a PDE model of          %
%%% evolutionary dynamics of a well-mixed population of aggressive        %
%%% breast cancer cells from in vitro data on MCF7-sh-WISP2 cell line     %
%%% and bootstrapping for uncertainty quantification                      %
%%% Copyright (C) 2023 C. Villa                                           %
%%%                                                                       %
%%% GNU General Public License:  <https://www.gnu.org/licenses/>          %
%%%                                                                       %
%%% See `main_Calibration.m or `main_UQ.m' files for more details.        %
%%%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Function to choose parameter domain bounds                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [theta_m, theta_M] = fun_ParDomain

    %%% Initial domain assumed for each parameter 
    gammag_dom = [0.001,3];
    gammal_dom = [0.001,3];
    alphaG_dom = [0.01,10];
    alphaL_dom = [0.01,10];
    d_dom = [1e-08,1e-05];
    kappaG_dom = [1e-07,1e-05];
    kappaL_dom = [1e-07,1e-05];
    etaL_dom = [1e-12,1e-04];
    beta_dom = [1e-04,1e-01];
    lambdaL_dom = [0,1];
    yL_dom = [0,15];
    yH_dom = [35,100];
    zeta_dom = [0,100];
    nh_dom = [0.9,4];
    mh_dom = [0.9,4];    
    
    %lambdaG_dom = [0,1];
    %Gm_dom = [0,25];
    
    
    theta_dom = [gammag_dom; gammal_dom; alphaG_dom; alphaL_dom; d_dom; ...
               kappaG_dom; kappaL_dom; etaL_dom; beta_dom; lambdaL_dom; ...
               yL_dom; yH_dom; zeta_dom; nh_dom; mh_dom];
    
    theta_m = theta_dom(:,1);
    theta_M = theta_dom(:,2);

end
