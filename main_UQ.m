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
%                                                                         %
% Code for uncertainty quantification using bootstrapping                 %
%  - Random data resampling, with replacement, without repetition         %
%  - MLE found with recursive algorithm (in domain shrunk from original   %
%    large one, closer to MLE obtained using mean data, to minimse cost)  %
%  - Post-processing: calculate bootstrap statistics (mean, standard      %
%    deviation, BIAS with OPS and empirical 95% confidence interval)      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load MCF7-sh-WISP2 data and all possible resampling 
[Ud1] = fun_ExpData('exp1');    
[Ud2] = fun_ExpData('exp2');  
% [perm_U] = fun_resample(Ud1,Ud2);
% save('ResampledData.mat','perm_U');
load('ResampledData.mat')

%% Calculate and store MLE for each resampled data set
[Udsd] = fun_ExpData('stdev');   % Load standard deviation (for likelihood)
[Udm] = fun_ExpData('mean');     % Load mean (for likelihood)
[y_m, y_M] = fun_ParDomain;      % Load parameter domain
Niter = 1;                       % Number of iterative shrinks of domain
%%% Initialise storing variables and files
file_name_temp = 'Bootstrapping/LatestMLE.mat'; 
yval_MLE_stored = [];
y_MLE_stored = [];
i1 = 1;
%%% Alt: load saved results
load('ResampledMLE_par.mat')
i1 = size(y_MLE_stored,1) +1
%% Iterate and get bootstrapping MLE
for i=i1:size(perm_U,1)
    tic
    [yval_MLE,y_MLE] = fun_iterMLE(y_m,y_M,Niter,perm_U(i,:),Udsd,file_name_temp);
    yval_MLE_stored = [yval_MLE_stored; yval_MLE];
    y_MLE_stored = [y_MLE_stored; y_MLE];
    save('ResampledMLE_val.mat','yval_MLE_stored');
    save('ResampledMLE_par.mat','y_MLE_stored');
    i
    toc
    
end

%% Post-processing
[y_m, y_M] = fun_ParDomain;      % Load parameter domain
load('ResampledMLE_par.mat')
y_MLE_stored = [y_MLE_stored(1:200,:)] ;

%%% Calculate summary statistics
Nb = size(y_MLE_stored,1)
Np = size(y_MLE_stored,2)
BTS_mean = sum(y_MLE_stored,1)./Nb;
BTS_var = sum((y_MLE_stored-BTS_mean).^2,1)./(Nb-1);
BTS_sd = sqrt(BTS_var);
CI = 0.95;

%% Plot distributions
p_row = 4;
p_col = 4;
bxln = 0.8;

par_lab = {'$\gamma_G$','$\gamma_L$','$\alpha_G$','$\alpha_L$','$d$','$\kappa_G$','$\kappa_L$','$\eta_L$','$\beta$','$\lambda_L$','$y_L$','$y_H$','$\zeta$','$m$','$c$'};

load('OPS_MAIN_0921.mat') % MAIN (both diffusion and drift)
BIAS = y-BTS_mean;
res = ones(1,Np);
j = Nb;

for i=1:Np

    %%% Confidence interval
    y_sorted = res(i)*sort(y_MLE_stored(1:j,i));
    CI_lower(i) = y_sorted(max(round(j*(1-CI)/2),1));
    CI_upper(i) = y_sorted(round(j*(CI+(1-CI)/2)));
    CI_lower(i)
    CI_upper(i)
    %%% Probability distribution
    [f, xi] = ksdensity(res(i)*y_MLE_stored(1:j,i));
    f = f(find(xi>=0):end); % Impose positivity constraint of parameters
    xi = xi(find(xi>=0):end); 
    f = f/trapz(xi, f);
    bxwdt = 0.05*max(f);
    bxpos = 1.25*max(f);
    %%% Plot
    subplot(p_row,p_col,i)
    % %%%% Option 1: Normalise histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % histogram(y_sorted, 'NumBins', 20, 'Normalization', 'pdf', 'FaceAlpha', 0.2);
    % hold on
    % plot(xi,f,'b', 'LineWidth', bxln)
    %%%% Option 2: Separate y scales %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yyaxis right;
    histogram(y_sorted, 'NumBins', 15, 'FaceAlpha', 0.2);
    rightYLabel = ylabel('Histogram PDF');
    rightYLabel.FontSize = 12;
    ylim([0, 50]);
    hold on
    yyaxis left;
    plot(xi,f,'b', 'LineWidth', bxln)
    leftYLabel = ylabel('Kernel Density Units');
    leftYLabel.FontSize = 12;
    hold on
    x_fill = [xi(xi >= CI_lower(i) & xi <= CI_upper(i)), CI_upper(i), CI_lower(i)];
    f_fill = [f(xi >= CI_lower(i) & xi <= CI_upper(i)), 0, 0];
    fill(x_fill, f_fill, 'b', 'FaceAlpha', 0.2);  % Transparent orange fill
    scatter(res(i)*BTS_mean(i), bxpos, 'ro', 'LineWidth', bxln, 'MarkerFaceColor', 'r');
    plot(res(i)*[BTS_mean(i)-BTS_sd(i),BTS_mean(i)+BTS_sd(i)], [bxpos,bxpos], '-r', 'LineWidth', 2*bxln);
    xline(res(i)*y(i),'g','LineWidth', 2)
    ylim([0,1.5*max(f)])
    title(par_lab(i), 'Interpreter','latex', 'FontSize', 18)
    axis square

end
lgd = legend({'KDE', '95% CI', 'Mean', 'Stdvn','MLE','Hist'}, 'Location', 'best', 'FontSize', 18,'Location','East');
lgd.Position(1) = 0.8;
subplot(p_row,p_col,11)


%% Simulate from Empirical 95% CI

%%% Sample parameter values from KDE
Ns = 200;     
par_s = zeros(Ns,Np);
for i=1:Np
    i
    samples = [];
    [f, xi] = ksdensity(y_MLE_stored(:,i));
    while numel(samples) < Ns
        p_s = randsample(xi, Ns-numel(samples), true, f);
        samples = [samples, p_s(p_s >= CI_lower(i) & p_s <= CI_upper(i))];
    end
    par_s(1:Ns,i) = samples;
end

%%% Simulate and save summary statistics to plot
[Udm] = fun_ExpData('mean');     
IC = [Udm(1),Udm(6),Udm(11),Udm(16)];
rhovs = [];
Gvs = [];
Lvs = [];
muvs = [];
muvr = [];
sig2vs = [];
for i=1:Ns
    i
    [tvs,rhos,Gs,Ls,mus,sig2s] = Extra_fun_SimPDE_toplot(IC,par_s(i,:));
    
    rhovs = [rhovs;rhos];
    Gvs = [Gvs;Gs];
    Lvs = [Lvs;Ls];
    muvs = [muvs;mus(1:length(tvs))];
    muvr = [muvr;mus(length(tvs)+1:end)];
    sig2vs = [sig2vs;sig2s];
end
tvr = tvs(end-size(muvr,2)+1:end);
