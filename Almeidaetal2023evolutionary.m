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
% %                                                                       %
% % Code to calibrate the partial differential equation model proposed    %
% % in Almeidaetal2022evolutionary (hereon referred to as "Al22Ev")       %
% % with the experimental data on MCF7-sh-WISP2 cells.                    %
% % Data on cell numbers, glucose and lactate concentrations, and mean    %
% % levels of MCT1 gene expression (fluorescence intensity from flow      %
% % cytometry experiments) is provided over the span of 5 days.           %
% % Model calibration relies on the function bayesopt (based on           %
% % bayesian optimisation) to maximise the likelyhood of a parameter      %
% % set drawn from a given bounded domain. The likelyhood is computed     %
% % assuming Gaussian measurement noise with zero mean, minimising the    %
% % distance between the normalised data and related summary statistics   %
% % obtained from simulating the nondimensional version of model given    %
% % coupled with Eq.(S18), Eq.(S19) and initial conditions (S22)-(S24).   %
% % For more details on the model nondimensionalisation, calibration      %
% % method, numerical scheme and biological problem see Al22Ev and        %
% % related Supplementary Material (also available on ArXiv).             %
% %                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                       %
%%% Almeidaetal2023evolutionary.m: calibration of a PDE model of        %
%%% evolutionary dynamics of a well-mixed population of aggressive        %
%%% breast cancer cells from in vitro data on MCF7-sh-WISP2 cell line     %
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
clc
clear all
close all
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 16)
set(0,'defaultlinelinewidth',2)

%% Data MCF7-sh-WISP2 (cf. Figure 1, Figure 2 of Al22Ex)

%%% Main experiment (1.5M, 1g/l) + Rescue
rhod = [1500,6875,6634,6003,7195]*1.0e+03;  % Cell number 
Gd = 1-[0.00,0.60,1.00,0.99,1.00];          % Glucose concentration
Ld = [0.00,8.35,11.35,10.48,9.66];          % Lactate concentration
mud = [15.57,25.72,26.20,34.00];            % Mean MCT1 expression
mud = [mud,19.28,11.92];                    % Rescue data

%%% Nondimensionalise and store
NDp = max(rhod); 
NDG = 1; 
NDL = max(Ld); 
ND = [NDp,NDG,NDL];
Ud = [rhod/NDp,Gd/NDG,Ld/NDL];

%% Calibration

%%% Domain assumed for each ND parameter (initial guess)
% Maximum proliferation rate via glycolysis
gammag = optimizableVariable('gammag',[0,10]);
% Maximum proliferation rate reusing lactate
gammal = optimizableVariable('gammal',[0,10]);
% Glucose concentration at half receptor occupancy
alphaG = optimizableVariable('alphaG',[0,10]);
% Lactate concentration at half receptor occupancy
alphaL = optimizableVariable('alphaL',[0,10]);
% Rate of cell death due to intracellular competition
d = optimizableVariable('d',[0,10]);
% Glucose consumption rate
kappaG = optimizableVariable('kappaG',[0,10]);
% Lactate production rate
kappaL = optimizableVariable('kappaL',[0,10]);
% Conversion factor for lactate consumption 
etaL = optimizableVariable('etaL',[0,10]);
% Minimum rate of spontaneous changes in MCT1 expression
beta = optimizableVariable('beta',[0,1e-01]);  
% Maximum rate of environment-induced increase in MCT1 expression
lambdaL = optimizableVariable('lambdaL',[0,1]); 
% Rate of environment-induced decrease in MCT1 expression
lambdaG = optimizableVariable('lambdaG',[0,1]); 
% MCT1 level corresponding to the maximum rate of proliferation via glycolysis
yL = optimizableVariable('yL',[0,15]);
% MCT1 level corresponding to the maximum rate of proliferation via lactate reuse
yM = optimizableVariable('yM',[35,100]);
% Lactate-dependency coefficient of the rate of spontaneous changes in MCT1 expression
zeta = optimizableVariable('zeta',[0,100]);      
% Hill coefficient for glucose ligand-receptor dynamics
nh = optimizableVariable('nh',[0.9,4]);
% Hill coefficient for lactate ligand-receptor dynamics
mh = optimizableVariable('mh',[0.9,4]);
% Threshold glucose concentration for lactate uptake G^*
Gm = optimizableVariable('Gm',[0,4.5]);

%%% Domain assumed for each parameter (progressively updated guess)
gammag = optimizableVariable('gammag',[1.8,2]);
gammal = optimizableVariable('gammal',[1.8,2]);
alphaG = optimizableVariable('alphaG',[0.5,0.7]);
alphaL = optimizableVariable('alphaL',[0.5,0.7]);
d = optimizableVariable('d',[0.8,1]);
kappaG = optimizableVariable('kappaG',[1.2,1.5]);
kappaL = optimizableVariable('kappaL',[1.2,1.5]);
etaL = optimizableVariable('etaL',[0.1,0.2]);
beta = optimizableVariable('beta',[0.0002,0.0008]); 
lambdaL = optimizableVariable('lambdaL',[0.05,0.1]); 
lambdaG = optimizableVariable('lambdaG',[0.2,0.4]); 
yL = optimizableVariable('yL',[0,5]);
yM = optimizableVariable('yH',[65,70]);
zeta = optimizableVariable('zeta',[5,8]);
nh = optimizableVariable('nh',[0.9,1.1]);
mh = optimizableVariable('mh',[0.9,1.1]);
Gm = optimizableVariable('Gm',[0.9,1.1]);

%%% Calibrate by minimising -(Likelyhood)
fun = @(x) Likelyhood(x.gammag,x.gammal,x.alphaG,x.alphaL,x.d,x.kappaG,x.kappaL,x.etaL,x.beta,x.lambdaL,x.lambdaG,x.yL,x.yH,x.zeta,x.nh,x.mh,x.Gm,Ud,mud);
results = bayesopt(fun,[gammag,gammal,alphaG,alphaL,d,kappaG,kappaL,etaL,beta,lambdaL,lambdaG,yL,yH,zeta,nh,mh,Gm],'IsObjectiveDeterministic',true,'MaxObjectiveEvaluations',100)

%%% Plot 
Simulation_plot(results.XAtMinObjective{:,1:17},Ud,mud,ND) 
y = results.XAtMinObjective{:,1:17}
yval = results.MinEstimatedObjective;

save('Optimal_1901_calib01.mat','yval','y')


%% ANNEXED FUNCTIONS - compute Likelyhood 

function mLD = Likelyhood(gammag,gammal,alphaG,alphaL,d,kappaG,kappaL,etaL,beta,lambdaL,lambdaG,ym,yM,zeta,hn,hm,Gm,Ud,mud)

    %%% Normalise Mean - [EQ. S2]
    mud = (mud-ym)/(yM-ym); 
    
    %%% Reproduce Rescue experiment if data is given as input
    Rescue = (length(mud)>4); 
    
    %%% Simulate and get summary statistics (given rho0 and G0)
    [Ui,mui] = Simulation(Ud(1),Ud(6),[gammag,gammal,alphaG,alphaL,d,kappaG,kappaL,etaL,beta,lambdaL,lambdaG,ym,yM,zeta,hn,hm,Gm],Rescue);

    %%% Compute likelyhood (short of 1/sqrt(2pi*LDvar) factor) - [EQ. S26]
    LD = 1;
    LDvar = 10; % Variance of Gaussian Likelyhood
    for k=1:length(mud)
        LD = LD*exp(-((mud(k)-mui(k))^2)/(2*LDvar));
    end
    for k=1:length(Ud)
        LD = LD*exp(-((Ud(k)-Ui(k))^2)/(2*LDvar));
    end 
    mLD = -LD; % To maximise likelyhood, minimise -LD
    
end

%% ANNEXED FUNCTIONS - simulate and return summary statistics

function [Ustore,mustore] = Simulation(rho0,G0,y,Rescue,NRescued)

    %%% Assign parameter values
    gammag = y(1);         % Maximum proliferation rate via glycolysis
    gammal = y(2);         % Maximum proliferation rate reusing lactate
    alphaG = y(3);         % Glucose concentration at half receptor occupancy
    alphaL = y(4);         % Lactate concentration at half receptor occupancy
    d = y(5);              % Rate of cell death due to intracellular competition
    kappaG = y(6);         % Glucose consumption rate
    kappaL = y(7);         % Lactate production rate
    etaL = y(8);           % Conversion factor for lactate consumption
    beta = y(9);           % Minimum rate of spontaneous changes in MCT1 expression
    lambdaL = y(10);       % Maximum rate of environment-induced increase in MCT1 expression
    lambdaG = y(11);       % Rate of environment-induced decrease in MCT1 expression
    yL = y(12);            % MCT1 level corresponding to the maximum rate of proliferation via glycolysis
    yH = y(13);            % MCT1 level corresponding to the maximum rate of proliferation via lactate reuse
    zeta = y(14);          % Lactate-dependency coefficient of the rate of spontaneous changes in MCT1 expression
    nh = y(15);            % Hill coefficient for glucose ligand dynamics
    mh = y(16);            % Hill coefficient for lactate ligand dynamics
    Gm = y(17);            % Threshold glucose concentration for lactate uptake G^*
    
    %%% Discretisation
    dt = 0.0001;
    t = 0:dt:5;
    Nt = length(t);
    xmi = -2;
    xm = 3;
    Nx = 1000; 
    x = linspace(xmi, xm, Nx);
    dx = abs(x(2)-x(1));

    %%% Initial conditions - [Eq. S21, S22, S23, S24]
    G = G0;                    % Initial glucose concentration
    L = 0;                     % Initial lactate concentration
    
    % % % % % % % Gaussian (with mean 15.57) 
    xi0 = (15.57-yL)/(yH-yL);  % Initial location (mean)        
    omega20 = 8/(yH-yL)^2;     % Initial squared scale (variance)
    alpha0 = 0;                % Initial SG parameter (no skewedness)   
    n = SKEWED_G(x,xi0,sqrt(omega20),alpha0,rho0);
    
    % % % % % % % Skewed Gaussian (with mean 15.57) % % % % %
    % xi0 = (11.01-ym)/(yM-ym);  % Initial location         %
    % omega20 = 35/(yM-ym)^2;    % Initial squared scale    %
    % alpha0 = 3.8;              % Initial SG parameter     %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    %%% Moments of the distribution - [EQ. S5_2]
    % Integral terms: left Riemann sum approximation
    rho = sum(n(1:end-1))*dx;                    % Cell number 
    mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;   % Mean trait

    %%% Matrix to store data
    Ustore = zeros(1,15);  
    Ustore(1) = rho0;
    Ustore(6) = G0;
    Ustore(11) = 0;
    mustore(1) = mu;
    tstore = [2,3,4,5];
    
    %%% p_L=l*U_L - [Eq. S8]
    g = gammag*(1-x.^2);
    %%% p_L=l*U_L - [Eq. S8]
    l = gammal*(1-(1-x).^2);
    lp = max(gammal*(1-(1-x).^2),0);

    %%% Iterate over time and compute approximate numerical solution
    for i=1:Nt
        
        %%% Heaviside step function H(G^*-G) - [Eq. S10]
        if G>Gm
            up = 0;
        else
            up = 1;
        end
        
        %%% Hill function (ligand dynamics) - [Eq. S9]
        Gup = G^nh/(alphaG^nh+G^nh);
        Lup = up*L^mh/(alphaL^mh+L^mh);
        
        %%% Fitness function R - [Eq. S6, S7, S8]
        R = g*Gup + l*Lup - d*rho;
        
        %%% Diffusion term + Phi - [Eq. S5_1, S15]
        % Second order central difference approximation
        Phi = beta*(1+zeta*Lup);
        diff = Phi*(n(3:Nx)+n(1:Nx-2)-2*n(2:Nx-1))/dx^2; 
        
        %%% Advection term + Psi - [Eq. Eq. S5_1, S16, S17]
        % First order upwind scheme
        Psi = lambdaL*Lup - (1-up)*lambdaG*max(0,mu);
        drift = (max(0,Psi)*(n(2:Nx-1)-n(1:Nx-2))+min(0,Psi)*(n(3:Nx)-n(2:Nx-1)))/dx;
      
        %%% Cell dynamics - [Eq. S5_1]
        nip1(2:Nx-1) = n(2:Nx-1) + dt*(R(2:Nx-1).*n(2:Nx-1) + diff - drift);

        %%% Zero-flux boundary conditions - [Sec. S2.4]
        nip1(1) = (Phi/(Phi+Psi*dx))*nip1(2);
        nip1(Nx) = (Phi/(Phi-Psi*dx))*nip1(Nx-1);

        %%% Glucose dynamics - [Eq. S18]
        G = G + dt*(-kappaG*Gup*rho);
        %%% Lactate dynamics - [Eq. S19]
        L = L + dt*(kappaL*Gup*rho - etaL*Lup*sum(lp(1:end-1).*n(1:end-1))*dx);
        
        %%% Update storing variables  
        n = nip1;
        rho = sum(n(1:end-1))*dx;
        mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;
        
        %%% Store at data timepoints
        if any(tstore==t(i))
            Ustore(t(i)) = rho;
            Ustore(t(i)+5) = G;
            Ustore(t(i)+10) = L;
            if t(i)>2 
                mustore = [mustore,(sum(n(1:end-1).*x(1:end-1))*dx)/rho];
            end
        end
        
        %%% Save phenotypic distribution at day 3 for Rescue experiment
        if Rescue & t(i)==3
            iR = i;
            nR = n;
        end
            
    end
    
    %%% Reproduce Rescue experiment
    if Rescue
        n = nR;
        rho = sum(n(1:end-1))*dx;
        mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;
        G = G0*4.5;
        L = 0;
        
        for i=(iR+1):Nt
        
            %%% Heaviside step function H(G^*-G) - [Eq. S10]
            if G>Gm
                up = 0;
            else
                up = 1;
            end

            %%% Hill function (ligand dynamics) - [Eq. S9]
            Gup = G^nh/(alphaG^nh+G^nh);
            Lup = up*L^mh/(alphaL^mh+L^mh);

            %%% Fitness function R - [Eq. S6, S7, S8]
            R = g*Gup + l*Lup - d*rho;

            %%% Diffusion term + Phi - [Eq. Eq. S5_1, S15]
            % Second order central difference approximation
            Phi = beta*(1+zeta*Lup);
            diff = Phi*(n(3:Nx)+n(1:Nx-2)-2*n(2:Nx-1))/dx^2; 

            %%% Advection term + Psi - [Eq. Eq. S5_1, S16, S17]
            % First order upwind scheme
            Psi = lambdaL*Lup - (1-up)*lambdaG*max(0,mu);
            drift = (max(0,Psi)*(n(2:Nx-1)-n(1:Nx-2))+min(0,Psi)*(n(3:Nx)-n(2:Nx-1)))/dx;

            %%% Cell dynamics - [Eq. S5_1]
            nip1(2:Nx-1) = n(2:Nx-1) + dt*(R(2:Nx-1).*n(2:Nx-1) + diff - drift);

            %%% Zero-flux boundary conditions - [Sec. S2.4]
            nip1(1) = (Phi/(Phi+Psi*dx))*nip1(2);
            nip1(Nx) = (Phi/(Phi-Psi*dx))*nip1(Nx-1);

            %%% Glucose dynamics - [Eq. S18]
            G = G + dt*(-kappaG*Gup*rho);
            %%% Lactate dynamics - [Eq. S19]
            L = L + dt*(kappaL*Gup*rho - etaL*Lup*sum(lp(1:end-1).*n(1:end-1))*dx);

            %%% Update storing variables  
            n = nip1;
            rho = sum(n(1:end-1))*dx;
            mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;
        
            %%% Store at data timepoints
            if t(i)==4 || t(i)==5
                mustore = [mustore,(sum(n.*x)*dx)/rho];
            end
        end
    end
    
end

%% ANNEXED FUNCTIONS - simulate and plot

function Simulation_plot(y,Ud,mud,ND)

    %%% Assign parameter values
    gammag = y(1);         % Maximum proliferation rate via glycolysis
    gammal = y(2);         % Maximum proliferation rate reusing lactate
    alphaG = y(3);         % Glucose concentration at half receptor occupancy
    alphaL = y(4);         % Lactate concentration at half receptor occupancy
    d = y(5);              % Rate of cell death due to intracellular competition
    kappaG = y(6);         % Glucose consumption rate
    kappaL = y(7);         % Lactate production rate
    etaL = y(8);           % Conversion factor for lactate consumption
    beta = y(9);           % Minimum rate of spontaneous changes in MCT1 expression
    lambdaL = y(10);       % Maximum rate of environment-induced increase in MCT1 expression
    lambdaG = y(11);       % Rate of environment-induced decrease in MCT1 expression
    yL = y(12);            % MCT1 level corresponding to the maximum rate of proliferation via glycolysis
    yH = y(13);            % MCT1 level corresponding to the maximum rate of proliferation via lactate reuse
    zeta = y(14);          % Lactate-dependency coefficient of the rate of spontaneous changes in MCT1 expression
    nh = y(15);            % Hill coefficient for glucose ligand dynamics
    mh = y(16);            % Hill coefficient for lactate ligand dynamics
    Gm = y(17);            % Threshold glucose concentration for lactate uptake G^*
    
    %%% Normalise Mean - [Eq. S2]
    mud = (mud-ym)/(yM-ym); 

    %%% Discretisation
    dt = 0.0001;
    tend = 5;
    t = 0:dt:tend;
    Nt = length(t);
    xmi = -2;
    xm = 3;
    Nx = 1000; 
    x = linspace(xmi, xm, Nx);
    dx = abs(x(2)-x(1));
    
    %%% Reproduce Rescue experiment if data is given as input
    Rescue = (length(mud)>4); 

    %%% Initial conditions - [Eq. S21, S22, S23, S24]
    G = Ud(6);                 % Initial glucose concentration
    L = 0;                     % Initial lactate concentration
    
    % Gaussian (with mean 15.57) 
    rho0 = Ud(1);              % Initial cell number
    xi0 = (15.57-ym)/(yM-ym);  % Initial location (mean)        
    omega20 = 8/(yM-ym)^2;     % Initial squared scale (variance)
    alpha0 = 0;                % Initial SG parameter (no skewedness)   
    n = SKEWED_G(x,xi0,sqrt(omega20),alpha0,rho0);
    
    %%% Moments of the distribution - [EQ. S5_2]
    % Integral terms: left Riemann sum approximation
    rho = sum(n(1:end-1))*dx;                    % Cell number 
    mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;   % Mean trait
    
    %%% Storing vectors 
    tstore = [1,2,3,4,5];
    nvs(1,:) = n;
    rhovs(1) = rho;
    Gvs(1) = G;
    Lvs(1) = L;
    tvs(1) = 0;
    muvs(1) = mu;
    sig2vs(1) = (sum(n(1:end-1).*(x(1:end-1).^2))*dx)/rho - mu^2;

    %%% p_L=l*U_L - [Eq. S8]
    g = gammag*(1-x.^2);
    %%% p_L=l*U_L - [Eq. S8]
    l = gammal*(1-(1-x).^2);
    lp = max(gammal*(1-(1-x).^2),0);

    %%% Iterate over time and compute approximate numerical solution
    for i=1:Nt
        
        %%% Heaviside step function H(G^*-G) - [Eq. S10]
        if G>Gm
            up = 0;
        else
            up = 1;
        end

        %%% Hill function (ligand dynamics) - [Eq. S9]
        Gup = G^nh/(alphaG^nh+G^nh);
        Lup = up*L^mh/(alphaL^mh+L^mh);

        %%% Fitness function R - [Eq. S6, S7, S8]
        R = g*Gup + l*Lup - d*rho;

        %%% Diffusion term + Phi - [Eq. S5_1, S15]
        % Second order central difference approximation
        Phi = beta*(1+zeta*Lup);
        diff = Phi*(n(3:Nx)+n(1:Nx-2)-2*n(2:Nx-1))/dx^2; 

        %%% Advection term + Psi - [Eq. S5_1, S16, S17]
        % First order upwind scheme
        Psi = lambdaL*Lup - (1-up)*lambdaG*max(0,mu);
        drift = (max(0,Psi)*(n(2:Nx-1)-n(1:Nx-2))+min(0,Psi)*(n(3:Nx)-n(2:Nx-1)))/dx;

        %%% Cell dynamics - [Eq. S5_1]
        nip1(2:Nx-1) = n(2:Nx-1) + dt*(R(2:Nx-1).*n(2:Nx-1) + diff - drift);

        %%% Zero-flux boundary conditions - [Sec. S2.4]
        nip1(1) = (Phi/(Phi+Psi*dx))*nip1(2);
        nip1(Nx) = (Phi/(Phi-Psi*dx))*nip1(Nx-1);

        %%% Glucose dynamics - [Eq. S18]
        G = G + dt*(-kappaG*Gup*rho);
        %%% Lactate dynamics - [Eq. S19]
        L = L + dt*(kappaL*Gup*rho - etaL*Lup*sum(lp(1:end-1).*n(1:end-1))*dx);

        %%% Update storing variables  
        n = nip1;
        rho = sum(n(1:end-1))*dx;
        mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;
        
        %%% Store evolution over time
        rhovs(i+1) = rho;
        Gvs(i+1) = G;
        Lvs(i+1) = L;
        muvs(i+1) = mu;
        sig2vs(i+1) = (sum(n(1:end-1).*(x(1:end-1).^2))*dx)/rho -(muvs(end))^2;
        tvs(i+1) = t(i);
        if any(tstore==t(i))
            nvs = [nvs;n];
        end
        
        if Rescue & t(i)==3
            iR = i;
            nR = n;
        end
            
    end

    %%% Reproduce Rescue experiment
    if Rescue
        n = nR; 
        rho = sum(n(1:end-1))*dx;
        mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;
        G = 4.5;
        L = 0;
        
        tvr = t(iR);
        muvr = mu;
        
        for i=(iR+1):Nt
        
            %%% Heaviside step function H(G^*-G) - [Eq. S10]
            if G>Gm
                up = 0;
            else
                up = 1;
            end

            %%% Hill function (ligand dynamics) - [Eq. S9]
            Gup = G^nh/(alphaG^nh+G^nh);
            Lup = up*L^mh/(alphaL^mh+L^mh);

            %%% Fitness function R - [Eq. S6, S7, S8]
            R = g*Gup + l*Lup - d*rho;

            %%% Diffusion term + Phi - [Eq. S5_1, S15]
            % Second order central difference approximation
            Phi = beta*(1+zeta*Lup);
            diff = Phi*(n(3:Nx)+n(1:Nx-2)-2*n(2:Nx-1))/dx^2; 

            %%% Advection term + Psi - [Eq. S5_1, S16, S17]
            % First order upwind scheme
            Psi = lambdaL*Lup - (1-up)*lambdaG*max(0,mu);
            drift = (max(0,Psi)*(n(2:Nx-1)-n(1:Nx-2))+min(0,Psi)*(n(3:Nx)-n(2:Nx-1)))/dx;

            %%% Cell dynamics - [Eq. S5_1]
            nip1(2:Nx-1) = n(2:Nx-1) + dt*(R(2:Nx-1).*n(2:Nx-1) + diff - drift);

            %%% Zero-flux boundary conditions - [Sec. S2.4]
            nip1(1) = (Phi/(Phi+Psi*dx))*nip1(2);
            nip1(Nx) = (Phi/(Phi-Psi*dx))*nip1(Nx-1);

            %%% Glucose dynamics - [Eq. S18]
            G = G + dt*(-kappaG*Gup*rho);
            %%% Lactate dynamics - [Eq. S19]
            L = L + dt*(kappaL*Gup*rho - etaL*Lup*sum(lp(1:end-1).*n(1:end-1))*dx);

            %%% Update storing variables  
            n = nip1;
            rho = sum(n(1:end-1))*dx;
            mu = (sum(n(1:end-1).*x(1:end-1))*dx)/rho;
        
            %%% Store at data timepoints
            tvr = [tvr,t(i)];
            muvr = [muvr,mu];
            if any(tstore==t(i))
                nvs = [nvs;n];
            end
        end
    end  

    %%% Plot solution
    %%% Rho
    set(gca,'LineWidth',3,'TickLength',[0.5 0.5]);
    subplot(2,3,1)
    plot(tvs,rhovs*ND(1))
    hold on
    scatter([0,1,2,3,4],Ud(1,1:5)*ND(1),60,'r','filled')
    legend('Model simulation','MCF7-sh-WISP2 data')
    grid on
    xlim([0,tend])
    ylim([0,max(8000000,max(rhovs)+10^4)])
    xticks([0,1,2,3,4,5])
    xlabel('$t$')
    axis square
    title(['$\rho(t)$'])
    %%% G
    subplot(2,3,2)
    plot(tvs,Gvs*ND(2))
    hold on
    scatter([0,1,2,3,4],Ud(1,6:10)*ND(2),60,'r','filled')
    grid on
    xlim([0,tend])
    legend('Model simulation','MCF7-sh-WISP2 data')
    xticks([0,1,2,3,4,5])
    xlabel('$t$')
    axis square
    title(['$G(t)$'])
    %%% L
    subplot(2,3,3)
    plot(tvs,Lvs*ND(3))
    hold on
    scatter([0,1,2,3,4],Ud(11:15)*ND(3),60,'r','filled')
    grid on
    xlim([0,tend])
    xticks([0,1,2,3,4,5])
    legend('Model simulation','MCF7-sh-WISP2 data')
    xlabel('$t$')
    axis square
    title(['$L(t)$'])
    %%% mu
    subplot(2,3,4)
    ym = 10^3*ym;
    yM = 10^3*yM;
    mud = mud*(yM-ym) + ym;  % Rescale back - [Eq. S4]
    muvs = muvs*(yM-ym) + ym;
    muvr = muvr*(yM-ym) + ym;
    plot(tvs(2:end),muvs(2:end))
    hold on
    scatter([0,3,4,5],mud(1:4),60,'r','filled')
    if Rescue
        hold on
        plot(tvr,muvr,'--','Color',[0 0.4470 0.7410])
        hold on
        scatter([4,5],mud(5:6),60,'r','filled','^')
    end
    xlim([0,tend])
    ylim([ym,yM])
    xticks([0,1,2,3,4,5])
    grid on
    xlabel('$t$')
    if Rescue
        legend('Model simulation','MCF7-sh-WISP2 data',"Model simulation" + newline  + "(rescue experiment)",'Rescue data')
    else
        legend('Model simulation','MCF7-sh-WISP2 data')
    end
    axis square
    title(['$\mu(t)$'])
    %%% sig2
    subplot(2,3,5)
    sig2vs = (sig2vs)*(yM-ym)^2;  % Rescale back - [Eq. S4]
    plot(tvs,sig2vs)
%     hold on
%     sig2est = [8,40,80,110]*10^6;    % Estimated variance
%     scatter([0,3,4,5],sig2est,60,'r','LineWidth',1.5)
    xlim([0,tvs(end)])
    grid on
    xticks([0,1,2,3,4,5])
%     legend('Model simulation','Est. from data')
    legend('Model simulation')
    xlabel('$t$')
    axis square
    title(['$\sigma^2(t)$'])
    
    %%% n
    x = (ym +x.*(yM-ym)); % Rescale back - [Eq. S2]
    posz = find(abs(x)==min(abs(x)))+1;
    x = x(posz:end);
    % Ignore rescaling factor in n because we are plotting it normalised
    nt = [nvs(1,posz:end)./max(nvs(1,:));nvs(4,posz:end)./max(nvs(4,:));nvs(5,posz:end)./max(nvs(5,:));nvs(6,posz:end)./max(nvs(6,:));nvs(7,posz:end)./max(nvs(7,:));nvs(8,posz:end)./max(nvs(8,:))];

    subplot(2,3,6)
    stepsize=0.6;
    bottom=0;
    %%% Rescue day 5
    semilogx(x,nt(6,:)+6*stepsize,'Color',[0.33, 0.42, 0.18]);
    hold on
    p = patch([x,fliplr(x)], [nt(6,:)+6*stepsize, repmat(bottom,size(x))], [0.33, 0.42, 0.18]);
    hold on
    data2 = 6*stepsize*ones(1,length(x))-0.01;
    area(x, data2, 'FaceColor', 'w','EdgeColor' , 'none');  
    %%% Rescue day 4
    semilogx(x,nt(5,:)+5*stepsize,'Color',[0.74, 0.2, 0.64]);
    p = patch([x,fliplr(x)], [nt(5,:)+5*stepsize, repmat(bottom,size(x))], [0.74, 0.2, 0.64]);
    hold on
    data2 = 5*stepsize*ones(1,length(x))-0.01;
    area(x, data2, 'FaceColor', 'w','EdgeColor' , 'none');
    %%% Main day 5
    semilogx(x,nt(4,:)+4*stepsize,'Color','y');
    p = patch([x,fliplr(x)], [nt(4,:)+4*stepsize, repmat(bottom,size(x))], 'y');
    hold on
    data2 = 4*stepsize*ones(1,length(x))-0.01;
    area(x, data2, 'FaceColor', 'w','EdgeColor' , 'none');
    %%% Main day 4
    semilogx(x,nt(3,:)+3*stepsize,'Color','gr');
    p = patch([x,fliplr(x)], [nt(3,:)+3*stepsize, repmat(bottom,size(x))], 'gr');
    hold on
    data2 = 3*stepsize*ones(1,length(x))-0.01;
    area(x, data2, 'FaceColor', 'w','EdgeColor' , 'none');
    %%% Main day 3
    semilogx(x,nt(2,:)+2*stepsize,'Color',[1.0, 0.55, 0.0]);
    p = patch([x,fliplr(x)], [nt(2,:)+2*stepsize, repmat(bottom,size(x))], [1.0, 0.55, 0.0]);
    hold on
    data2 = 2*stepsize*ones(1,length(x))-0.01;
    area(x, data2, 'FaceColor', 'w','EdgeColor' , 'none');
    %%% Main day 0
    semilogx(x,nt(1,:)+1*stepsize,'Color','b');
    p = patch([x,fliplr(x)], [nt(1,:)+1*stepsize, repmat(bottom,size(x))], 'b');
    hold on
    data2 = 1*stepsize*ones(1,length(x))-0.01;
    area(x, data2, 'FaceColor', 'w','EdgeColor' , 'none');

    xlim([10,1.68e+05])
    ylim([0 8*stepsize])
    set(gca, 'Layer', 'top')
    yticks([1 2 3 4 5 6]*stepsize)
    set(gca,'TickLabelInterpreter','latex')
    yticklabels({'$t=0$','$t=3$','$t=4$','$t=5$','$t=4$ R','$t=5$ R'})
    xticks([10^0,10^1,10^2,10^3,10^4,10^5,10^6])
%     xline(10^4,'-','Linewidth',1)
%     xline(10^5,'-','Linewidth',1) 
    xlabel('$y$')
    title('$n(t,y)$')
    hold off;
    axis square
    drawnow
    
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
