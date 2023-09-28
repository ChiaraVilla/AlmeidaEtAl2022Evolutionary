%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% "Experimentally-informed mathematical modelling of the evolutionary %%%
%%%   dynamics of glucose-deprived aggressive breast cancer cells"      %%%
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
% Function to simulate the PDE model developped in the above-mentioned    %
% manuscript and return the summary statistics to compare with the data   %
% for a give (input) parameter set                                        %
%                                                                         %
% INPUT:                                                                  %
%   - Initial condition (IC): cell number, glucose and lactate conc, mu   %
%   - Parameter set (y)                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Us] = fun_SimPDE(IC,y)

    %%% Assign parameter values
    gammag = y(1);         % Maximum proliferation rate via glycolysis
    gammal = y(2);         % Maximum proliferation rate reusing lactate
    alphaG = y(3);         % Glucose concentration at half receptor occupancy
    alphaL = y(4);         % Lactate concentration at half receptor occupancy
    d = y(5);              % Rate of cell death due to competition for space
    kappaG = y(6);         % Glucose consumption rate
    kappaL = y(7);         % Lactate production rate
    etaL = y(8);           % Conversion factor for lactate consumption
    beta = y(9);           % Minimum rate of phenotypic variation
    lambdaL = y(10);       % Maximum rate of lactate-driven MCT1 expression increase
    ym = y(11);            % Fluorescence intesity corresponding to x=0
    yM = y(12);            % Fluorescence intesity corresponding to x=1
    zeta = y(13);          % Lactate-induced phenotypic variation coefficient
    nh = y(14);            % Hill coefficient for glucose ligand dynamics
    mh = y(15);            % Hill coefficient for lactate ligand dynamics
    
    %%% Discretisation
    dt = 1e-05;
    t = 0:dt:5;
    Nt = length(t);
    xmi = -2;
    xm = 3;
    Nx = 1000; 
    x = linspace(xmi, xm, Nx);
    dx = abs(x(2)-x(1));

    %%% Initial conditions - [Eq. A.16, A.17, A.18]
    rho = IC(1);              % Initial cell number
    G = IC(2);                % Initial glucose concentration
    L = IC(3);                % Initial lactate concentration
    mu = IC(4);               % Initial mean MCT1 expression
    
    % % % % % % % Gaussian (with mean 15.57) 
    xi0 = (mu-ym)/(yM-ym);     % Initial location (mean)        
    omega20 = 8/(yM-ym)^2;     % Initial squared scale (variance)
    alpha0 = 0;                % Initial SG parameter (no skewedness)   
    n = SKEWED_G(x,xi0,sqrt(omega20),alpha0,rho);
    
    % % % % % % % Skewed Gaussian (with mean 15.57) % % % % %
    % xi0 = (11.01-ym)/(yM-ym);  % Initial location         %
    % omega20 = 35/(yM-ym)^2;    % Initial squared scale    %
    % alpha0 = 3.8;              % Initial SG parameter     %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    %%% Moments of the distribution - [EQ. A.2]
    % Integral terms: left Riemann sum approximation
    rho = sum(n(1:end-1))*dx;                    % Cell number 
    mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;   % Mean trait

    %%% Matrix to store data
    rhostore = rho;
    Gstore = G;
    Lstore = L;
    mustore = mu;
    tstore_main = [1,2,3,4];
    tstore_mu = [3,4,5];

    %%% p_L=l*U_L - [Eq. A.6]
    g = gammag*(1-x.^2);
    %%% p_L=l*U_L - [Eq. A.7]
    l = gammal*(1-(1-x).^2);
    lp = max(gammal*(1-(1-x).^2),0);

    %%% Iterate over time and compute approximate numerical solution
    for i=1:Nt
        
        %%% Hill function (ligand dynamics) - [Eq. A.8]
        Gup = G^nh/(alphaG^nh+G^nh);
        Lup = L^mh/(alphaL^mh+L^mh);
        
        %%% Fitness function R - [Eq. A.4, A.5, A.6, A.7]
        R = g*Gup + l*Lup - d*rho;
        
        %%% Diffusion term + Phi - [Eq. A.11]
        % Second order central difference approximation
        Phi = beta*(1+zeta*Lup);
        diff = Phi*(n(3:Nx)+n(1:Nx-2)-2*n(2:Nx-1))/dx^2; 
        
        %%% Advection term + Psi - [Eq. A.12, A.13]
        % First order upwind scheme
        Psi = lambdaL*Lup;
        drift = (max(0,Psi)*(n(2:Nx-1)-n(1:Nx-2))+min(0,Psi)*(n(3:Nx)-n(2:Nx-1)))/dx;
      
        %%% Cell dynamics - [Eq. A.3]
        nip1(2:Nx-1) = n(2:Nx-1) + dt*(R(2:Nx-1).*n(2:Nx-1) + diff - drift);

        %%% Zero-flux boundary conditions - [Eq. S3]
        nip1(1) = (Phi/(Phi+Psi*dx))*nip1(2);
        nip1(Nx) = (Phi/(Phi-Psi*dx))*nip1(Nx-1);

        %%% Glucose dynamics - [Eq. A.14]
        G = G + dt*(-kappaG*Gup*rho);
        if G<0
            G=0;
        end
        %%% Lactate dynamics - [Eq. A.15]
        L = L + dt*(kappaL*Gup*rho - etaL*Lup*sum(lp(1:end-1).*n(1:end-1))*dx);
        
        %%% Update storing variables  
        n = nip1;
        rho = sum(n(1:end-1))*dx;
        mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;
        
        %%% Store at data timepoints
        if any(tstore_main==t(i))
            rhostore = [rhostore,rho];
            Gstore = [Gstore,G];
            Lstore = [Lstore,L];
        end
        if any(tstore_mu==t(i))
            mustore = [mustore,mu*(yM-ym) + ym];
        end
        
    end

    Us = [rhostore,Gstore,Lstore,mustore];
    
end

%% ANNEXED FUNCTIONS - construct Weighted Skewed-Gaussian distribution

function n = SKEWED_G(x,eps,w,alp,rho)
    n = rho*(2/w)*(psi_G(x,eps,w)).*(PSI_cdf(x,eps,w,alp));
end
function psi = psi_G(x,eps,w)
    z = (x-eps)/w;
    psi = (1/sqrt(2*pi))*exp(-(z.^2)/2);
end
function PSI = PSI_cdf(x,eps,w,alp)
    z = alp*((x-eps)/w);
    PSI = 0.5*(1+erf(z/sqrt(2)));
end
