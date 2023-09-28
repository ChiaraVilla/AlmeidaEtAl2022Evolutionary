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
% Store experimental data for model calibration and uncertainty           %
% quantification: sample experiment 1 ('_1'), sample experiment 2 ('_2'), %
% mean value ('_m') and standard deviation ('_sd')                        %
%                                                                         %
% Input data_type required, options: 'exp1', 'exp2', 'mean', 'stdev'      %
% Data deturned in vector containting [rho,G,L,mu]                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data] = fun_ExpData(data_type)

%%% MCF7-sh-WISP2 cell number (1.5M, 1g/l glucose)
rhod_1 = [1500,6875,6634,6003,7195]*1e+03;
rhod_2 = [1500,5800,5600,5000,6200]*1e+03;
rhod_m = [1500,6337.5,6117,5501.5,6697.5]*1e+03;
rhod_sd = [0,760139.79,731148.41,709228.10,703571.25];

%%% Glucose concentration (1.5M, 1g/l glucose - in mM)
Gd_1 = [5.46,0,0,0,0];
Gd_2 = [5.58,1.11,0.02,0.06,0.01];
Gd_m = [5.52,0.555,0.01,0.03,0.005];
Gd_sd = [0.084852814,0.784888527,0.014142136,0.042426407,0.007071068];

%%% Lactate concentration (1.5M, 1g/l glucose - in mM)
Ld_1 = [1.67,9.26,7.94,6.58,5.51];
Ld_2 = [1.67,8.35,11.35,10.48,9.66];
Ld_m = [1.67,8.805,9.645,8.53,7.585];
Ld_sd = [0,0.643467171,2.411234124,2.757716447,2.934493142];

%%% MCT1 expression (1.5M, 1g/l glucose + Rescue 4.5g/l glucose)
%%% Glucose deprivation experiment days 0,3,4,5, Rescue days 4,5
mud_1 = [15.52,25.66,25.98,33.9,19.27,11.8];
mud_2 = [15.61,25.78,26.41,34.05,19.29,12.04];
mud_m = [15.57,25.72,26.20,33.98,19.28,11.92];
mud_sd = [0.06363961,0.084852814,0.304055916,0.106066017,0.014142136,0.169705627];

%%% Select data to return
switch data_type
    case 'exp1'
        data = [rhod_1,Gd_1,Ld_1,mud_1(1:4)];
    case 'exp2'
        data = [rhod_2,Gd_2,Ld_2,mud_2(1:4)];
    case 'mean'
        data = [rhod_m,Gd_m,Ld_m,mud_m(1:4)];
    case 'stdev'
        data = [rhod_sd,Gd_sd,Ld_sd,mud_sd(1:4)];
    otherwise
        error('Specify data type: exp1, exp2, mean, stdev')
end
