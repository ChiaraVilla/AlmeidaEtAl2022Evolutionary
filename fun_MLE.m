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
%  Function to calculate MLE given parameter domain bounds in input and   %
%  the maximum number of function evaluations using baysopt               %
%  Inputs:                                                                %
%  - y_m: lower bound of parameter space domain                           %
%  - y_M: upper bound of parameter space domain                           %
%  - Neval: maximum number of baysopt function evaluations                %
%  - Udm: data points to fit (required for Likelihood function)           %
%  - Udsd: variance of each data point (for Likelihood function)          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [yval,y_MLE] = fun_MLE(y_m,y_M,NEval,Udm,Udsd)
    
    % Maximum proliferation rate via glycolysis
    gammag = optimizableVariable('gammag',[y_m(1),y_M(1)]);
    % Maximum proliferation rate reusing lactate
    gammal = optimizableVariable('gammal',[y_m(2),y_M(2)]);
    % Glucose concentration at half receptor occupancy
    alphaG = optimizableVariable('alphaG',[y_m(3),y_M(3)]);
    % Lactate concentration at half receptor occupancy
    alphaL = optimizableVariable('alphaL',[y_m(4),y_M(4)]);
    % Rate of cell death due to intracellular competition
    d = optimizableVariable('d',[y_m(5),y_M(5)]);
    % Glucose consumption rate
    kappaG = optimizableVariable('kappaG',[y_m(6),y_M(6)]);
    % Lactate production rate
    kappaL = optimizableVariable('kappaL',[y_m(7),y_M(7)]);
    % Conversion factor for lactate consumption 
    etaL = optimizableVariable('etaL',[y_m(8),y_M(8)]);
    % Minimum rate of spontaneous changes in MCT1 expression
    beta = optimizableVariable('beta',[y_m(9),y_M(9)]);  
    % Maximum rate of environment-induced increase in MCT1 expression
    lambdaL = optimizableVariable('lambdaL',[y_m(10),y_M(10)]);  
    % MCT1 level corresponding to the maximum rate of proliferation via glycolysis
    yL = optimizableVariable('yL',[y_m(11),y_M(11)]);
    % MCT1 level corresponding to the maximum rate of proliferation via lactate reuse
    yH = optimizableVariable('yH',[y_m(12),y_M(12)]);
    % Lactate-dependency coefficient of the rate of spontaneous changes in MCT1 expression
    zeta = optimizableVariable('zeta',[y_m(13),y_M(13)]);      
    % Hill coefficient for glucose ligand-receptor dynamics
    nh = optimizableVariable('nh',[y_m(14),y_M(14)]);
    % Hill coefficient for lactate ligand-receptor dynamics
    mh = optimizableVariable('mh',[y_m(15),y_M(15)]);
    
    %%% Maximise Likelihood (minimise - Likelihood)
    fun = @(x) fun_mLikelyhood(x.gammag,x.gammal,x.alphaG,x.alphaL,x.d,x.kappaG,x.kappaL,x.etaL,x.beta,x.lambdaL,x.yL,x.yH,x.zeta,x.nh,x.mh,Udm,Udsd);
    results = bayesopt(fun,[gammag,gammal,alphaG,alphaL,d,kappaG,kappaL,etaL,beta,lambdaL,yL,yH,zeta,nh,mh],'IsObjectiveDeterministic',true,'MaxObjectiveEvaluations',NEval,'PlotFcn',[])
    y_MLE = results.XAtMinObjective{:,1:15}
    yval = results.MinEstimatedObjective;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Function to calculate - ("minus") log(Likelihood) to find MLE          %
%  Inputs:                                                                %
%  - parameter values                                                     %
%  - Udm: data points to fit                                              %
%  - Udsd: variance of each data point                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mLD = fun_mLikelyhood(gammag,gammal,alphaG,alphaL,d,kappaG,kappaL,etaL,beta,lambdaL,ym,yM,zeta,hn,hm,Udm,Udsd)

    %%% Initial conditions (rho, G, L, mu)
    IC = [Udm(1),Udm(6),Udm(11),Udm(16)];
 
    %%% Simulate and get summary statistics 
    y = [gammag,gammal,alphaG,alphaL,d,kappaG,kappaL,etaL,beta,lambdaL,ym,yM,zeta,hn,hm];
    [Ui] = fun_SimPDE(IC,y);

    %%% Compute likelyhood (ignore ICs) 
    Udm = [Udm(2:5),Udm(7:10),Udm(12:15),Udm(17:end)];
    Udsd = [Udsd(2:5),Udsd(7:10),Udsd(12:15),Udsd(17:end)];
    Ui = [Ui(2:5),Ui(7:10),Ui(12:15),Ui(17:end)];
    LD = 0;
    for k=1:length(Ui)
        %%% Weighted sum of squared residual
        % % (a) Set Gaussian error variance proportional to data value
        % LD = LD + (((Udm(k)-Ui(k))/Udm(k))^2)/2;
        % % (b) Set Gaussian error variance proportional to data variance
        LD = LD + (((Udm(k)-Ui(k))/Udsd(k))^2)/2;
    end
    mLD = LD;
    
end



