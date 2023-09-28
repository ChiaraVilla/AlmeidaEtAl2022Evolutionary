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
%%% This program is free software: you can redistribute it and/or modify  %
%%% it under the terms of the GNU General Public License as published by  %
%%% the Free Software Foundation, either version 3 of the License, or     %
%%% (at your option) any later version.                                   %
%%%                                                                       %
%%% This program is distributed in the hope that it will be useful,       %
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%%% GNU General Public License for more details.                          %
%%%                                                                       %
%%% You should have received a copy of the GNU General Public License     %
%%% along with this program.  If not, see <https://www.gnu.org/licenses/>.%
%%%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                                                       %
% % Code to calibrate the partial differential equation model proposed    %
% % in Almeidaetal2022evolutionary (hereon referred to as "Al22Ev")       %
% % with the experimental data on MCF7-sh-WISP2 cells.                    %
% % Data on cell numbers, glucose and lactate concentrations, and mean    %
% % levels of MCT1 gene expression (fluorescence intensity from flow      %
% % cytometry experiments) is provided over the span of 5 days.           %
% % Model calibration relies on the function bayesopt (based on           %
% % bayesian optimisation) to maximise the likelyhood of a parameter      %
% % set drawn from a given bounded domain. In particular, we minimise     %
% % the weighted sum of squared residuals of the summary statistics       %
% % obtained from simulating the nondimensional version of model given    %
% % by Eq.(S5) and definitions in Eq.(S6)-(S10) and (S15)-(S17),          %
% % coupled with Eq.(S18), Eq.(S19) and initial conditions (S22)-(S24),   %
% % capared witht he mean values of the summary statistics recorded       %
% % during replicate experiments. The weights are the inverse variance    %
% % calculated from experimental data to account for heteroschedasticity. %
% % For more details on the model nondimensionalisation, calibration      %
% % method, numerical scheme and biological problem see Al22Ev and        %
% % related Supplementary Material (also available on ArXiv).             %
% %                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Udsd] = fun_ExpData('stdev');   % Load standard deviation (for likelihood)
[Udm] = fun_ExpData('mean');     % Load mean (for likelihood)
[y_m, y_M] = fun_ParDomain;      % Load parameter domain
Niter = 1;                       % Number of iterative shrinks of domain

tic
[yval_MLE,y_MLE] = fun_iterMLE(y_m,y_M,Niter,Udm,Udsd,'OPS_main_0924.mat');
toc
