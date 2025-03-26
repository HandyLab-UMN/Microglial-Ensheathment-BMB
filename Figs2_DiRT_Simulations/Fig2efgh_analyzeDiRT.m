%%
% Analyze and plot the official results of the DiRT simulations (reproduces
% Fig 2e-g)
%%
clear; close all; clc;
restoredefaultpath;
folder = fileparts(which('Fig2efgh_analyzeDiRT.m')); 
addpath(genpath(folder));
rmpath(folder)
%% Load the averaged data for the non-symmetric case
load('./DiRT_Sim_Data/DiRT_Results_nonsym_N50_tau0.10_ave.mat')

%%
alphaFn = @(b,t) b(1)/b(2)^2*t.*exp(-t/b(2));
% First term estimated for phi = -1 and then fixed at 55.7345
modelFn3 = @(b,t) 55.7345*exp(-t/b(1));

numPhi = 1:size(phi_vec,2);

JMax = zeros(length(numPhi),1);
tau_s = zeros(length(numPhi),1);
rmse = zeros(length(numPhi),1);
rsquared = zeros(length(numPhi),1);

Jexp = zeros(length(numPhi),1);
tau_sExp = zeros(length(numPhi),1);
rmseExp = zeros(length(numPhi),1);
rsquaredExp = zeros(length(numPhi),1);

tic;
parfor ii = 1:length(numPhi)

    % Initial guesses [Jmax, tau_s]
    initialGuess = [max(DiRT_results_ave(:,2,ii))*exp(1)*0.1 0.1];

    % Fit the alpha-fn model
    alphaMDL = fitnlm(DiRT_results_ave(:,1,ii),DiRT_results_ave(:,2,ii),...
        alphaFn,initialGuess);


    expMDL = fitnlm(DiRT_results_ave(:,1,ii),DiRT_results_ave(:,2,ii),...
        modelFn3,initialGuess(2));

    % Save the values;
    JMax(ii) = alphaMDL.Coefficients.Estimate(1);
    tau_s(ii,1) = alphaMDL.Coefficients.Estimate(2);
    
    rmse(ii,:) = alphaMDL.RMSE;
    rsquared(ii,:) = alphaMDL.Rsquared.Ordinary;

    tau_sExp(ii) = expMDL.Coefficients.Estimate(1);
    rmseExp(ii,:) = expMDL.RMSE;
    rsquaredExp(ii,:) = expMDL.Rsquared.Ordinary;
end
toc;

%% Fit the piecewise linear functions

sympref('HeavisideAtOrigin', 0);

% Fit the piecewise linear function for J
pwlFn =  @(b,x) b(1)+(b(2).*(x-b(3))).*heaviside(x-b(3));
initialGuess = [JMax(1) -JMax(1) 0];
pwlMDL_J = fitnlm([phi_vec],[JMax],pwlFn,initialGuess);

% Fit the piecewse linear function for tau
notchLoc = pwlMDL_J.Coefficients.Estimate(3);
initialGuess = [tau_s(1) -tau_s(1)];
pwlFn_v2 =  @(b,x) b(1)+(b(2).*(x-notchLoc)).*heaviside(x-notchLoc);
pwlMDL_tau = fitnlm(phi_vec,tau_s,pwlFn_v2,initialGuess);

%% Plot the results

figure(1); clf; 
colorScheme = get(gca,'colororder');
exPlot1 = 1;
exPlot2 = 176;
tVec = [0:0.0001:1];

% Compare the simulations with the fits
subplot(1,5,1); hold on;
plot(DiRT_results_ave(:,1,exPlot1),DiRT_results_ave(:,2,exPlot1),...
    'color','k','LineWidth',1.5)
plot(tVec,alphaFn([JMax(exPlot1),tau_s(exPlot1)],tVec),...
    '--','LineWidth',1.5,'color','k')
plot(tVec,modelFn3([tau_sExp(exPlot1)],tVec),...
    '-.k','LineWidth',1.5)
xlim([0 0.75]);
ylim([0 60])
set(gca,'fontsize',16)
xlabel('Time (a.u.)')
ylabel('# Activated Receptors')
legend('Sim',sprintf('\\alpha-fn'),'Exp. fit')
title('\phi = -1')
yticks(0:10:60)

subplot(1,5,2); hold on;
plot(DiRT_results_ave(:,1,exPlot2),DiRT_results_ave(:,2,exPlot2),...
    'LineWidth',1.5,'color',colorScheme(1,:))
plot(tVec,alphaFn([JMax(exPlot2),tau_s(exPlot2)],tVec),...
    '--','LineWidth',1.5,'color',colorScheme(1,:))
plot(tVec,modelFn3([tau_sExp(exPlot2)],tVec),...
    '-.k','LineWidth',1.5,'color',colorScheme(1,:))
xlim([0 0.75]);
ylim([0 60])
title('\phi = 0.75')
set(gca,'fontsize',16)
xlabel('Time (a.u.)')
ylabel('# Activated Receptors')
yticks(0:10:60)

% Plot the R^2 values
subplot(1,5,3); hold on;
plot(phi_vec,rsquared,'k--','markersize',16,'LineWidth',1.5)
plot(phi_vec,rsquaredExp,'k-.','markersize',7,...
    'MarkerFaceColor',[0.8500    0.3250    0.0980],'LineWidth',1.5)
xlabel('\phi')
ylabel('R^2')
set(gca,'fontsize',16)
legend('\alpha function','Exp function')


% Plot all values \phi
subplot(1,5,4); hold on;
plot([phi_vec],[JMax],'.','markersize',15,'linewidth',1.5)
plot([phi_vec],pwlMDL_J.predict([phi_vec']),'k-','linewidth',1.5)
xlabel('Protrusion depth')
ylabel('J (synaptic strength)')
set(gca,'fontsize',16)
ylim([0 18])

subplot(1,5,5); hold on;
plot(phi_vec,tau_s,'.','markersize',15,'linewidth',1.5)
plot(phi_vec,pwlMDL_tau.predict(phi_vec'),'k-','linewidth',1.5)
set(gca,'fontsize',16)
xlabel('Protrusion depth')
ylabel('\tau_s (synaptic time constant)')
ylim([0.0 0.14])

fprintf('------\n');
fprintf('J Rsquared: %.4f\n',pwlMDL_J.Rsquared.ordinary)
fprintf('Notch location: %.4f\n', pwlMDL_J.Coefficients.Estimate(3))
fprintf('tau Rsquared: %.4f\n',pwlMDL_tau.Rsquared.ordinary)
fprintf('------\n');