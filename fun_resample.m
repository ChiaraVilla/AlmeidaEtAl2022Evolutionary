%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% "Experimentally-informed mathematical modelling of the evolutionary %%%
%%%   dynamics of glucose-deprived aggressive breast cancer cells"      %%%
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
% Function to create all possible resampling with replacement and without % 
% repetition fro bootstrapping algorithm (only for 2 data samples).       %
% Input: samples U1, U2.                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [perm_U] = fun_resample(U1,U2)

    N = length(U1); 
    
    perm_v = zeros(1,N);
    perm_U = U1;
    
    while size(perm_v, 1) < N^2
        rand_v = randi([0, 1], 1, N);
    
        if any(ismember(perm_v, rand_v, 'rows'))
            continue;  
        else
            perm_v = [perm_v; rand_v];
            perm_U = [perm_U; U1.*(rand_v==0) + U2.*(rand_v==1)];
        end
    end


end
