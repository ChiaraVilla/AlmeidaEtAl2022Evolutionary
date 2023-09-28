%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%       "Evolutionary dynamics of glucose-deprived cancer cells:      %%%
%%%     insights from experimentally-informed mathematical modelling"   %%%
%%%                                                                     %%%
%%%     L. Almeida, J. Denis, N. Ferrand, T. Lorenzi, M. Sabbah,        %%%
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
%  Function to calculate MLE given intial parameter domain bounds in      %
%  input and iteratively shrink the domain to improve the MLE estimate    %
%  Inputs:                                                                %
%  - y_m: initial lower bound of parameter space domain                   %
%  - y_M: initial upper bound of parameter space domain                   %
%  - Niter: maximum number of domain-shrinking iterations                 %
%  - Udm: data points to fit (required for Likelihood function)           %
%  - Udsd: variance of each data point (for Likelihood function)          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [yval_MLE,y_MLE] = fun_iterMLE(y_m,y_M,Niter,Udm,Udsd,file_name)

    %%% Start with many function evaluation points to estimate MLE
    Neval = [500,400,300,200,100,repelem(50,Niter-5)];
    yval_MLE = 5000;

    %%% Shrinking domain fraction
    sh_df = 2/3;

    %%% Desired shrinking factor
    sh_f = repelem(1,Niter);
    
    for i=1:Niter
        
        %%% Calculate MLE in given parameter domain
        [yval,y] = fun_MLE(y_m,y_M,Neval(i),Udm,Udsd);
        
        %%% Update MLE
        if yval<yval_MLE
            yval_MLE = yval
            y_MLE = y;
            save(file_name,'yval_MLE','y_MLE','y_m','y_M','Udm','Udsd')
        end
    
        %%% Shrink parameter domain 
        y_m = y_m + sh_df*(y_MLE'-y_m);
        y_M = y_M - sh_df*(y_M-y_MLE');
        %%% Logaritmic shrinking of parameter domain (if large domain)
        % y_m = max(y_m,(y_MLE').*exp(-sh_f(i)*(log(y_MLE')-log(y_m))));
        % y_M = min(y_M,(y_MLE').*exp(sh_f(i)*(log(y_M)-log(y_MLE'))));
    
    end

end