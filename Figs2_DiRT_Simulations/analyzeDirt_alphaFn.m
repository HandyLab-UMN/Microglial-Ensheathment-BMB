%%
% Analyze and plot the official results of the DiRT simulations (reproduces
% Fig 1)
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('analyzeDirt_alphaFn.m')); 
addpath(genpath(folder));
rmpath(folder)
%% Run DiRT_Simulation_Demo first
load('./DiRT_Sim_Data/DiRT_Results_N50_tau0.10_v2.mat')


%% Ave the results and get curve statistics

[DiRT_results_ave,alphaFn, JMax,tau_s,rmse,rsquared,...
    pwlMDL_J, pwlMDL_tau]=DiRT_Analysis_Fn(phi_vec,DiRT_results);
%%

figure(1); clf; 
subplot(1,2,1); hold on;
plot(phi_vec,JMax,'o','markersize',8,'linewidth',1.5)
plot(phi_vec,pwlMDL_J.predict(phi_vec'),'k-','linewidth',1.5)
xlabel('Protrusion depth')
ylabel('J (synaptic strength)')
set(gca,'fontsize',16)
ylim([0 18])

subplot(1,2,2); hold on;
plot(phi_vec,tau_s,'o','markersize',8,'linewidth',1.5)
plot(phi_vec,pwlMDL_tau.predict(phi_vec'),'k-','linewidth',1.5)
set(gca,'fontsize',16)
xlabel('Protrusion depth')
ylabel('\tau_s (synaptic time constant)')
ylim([0.02 0.14])

% plot([0:.01:2],pwlMDL_tau.predict([0:.01:2]'),'k--','linewidth',1.5)
fprintf('------\n');
fprintf('J Rsquared: %.4f\n',pwlMDL_J.Rsquared.ordinary)
fprintf('Notch location: %.4f\n', pwlMDL_J.Coefficients.Estimate(3))
fprintf('tau Rsquared: %.4f\n',pwlMDL_tau.Rsquared.ordinary)
fprintf('------\n');



%%
colorScheme = get(gca,'colororder');
exPlot1 = 1;
exPlot2 = find(phi_vec>=0.75,1);
tVec = [0:0.0001:1];
figure(2); clf; hold on;

plot(DiRT_results_ave(:,1,exPlot1),DiRT_results_ave(:,2,exPlot1),'color','k','LineWidth',1.5)
plot(tVec,alphaFn([JMax(exPlot1),tau_s(exPlot1)],tVec),'--','color','k','LineWidth',1.5)

plot(DiRT_results_ave(:,1,exPlot2),DiRT_results_ave(:,2,exPlot2),...
    'LineWidth',1.5,'color',colorScheme(1,:))
plot(tVec,alphaFn([JMax(exPlot2),tau_s(exPlot2)],tVec),...
    '--','LineWidth',1.5,'color',colorScheme(1,:))
set(gca,'fontsize',16)
xlim([0 1])
xlabel('Time (a.u.)')
ylabel('# Act. Rec.')


figure(3); clf; 
subplot(1,2,1); hold on;
plot(DiRT_results_ave(:,1,exPlot1),DiRT_results_ave(:,2,exPlot1),'color','k','LineWidth',1.5)
plot(DiRT_results_ave(:,1,exPlot2),DiRT_results_ave(:,2,exPlot2),...
    'LineWidth',1.5,'color',colorScheme(1,:))
xlim([0 0.75]);
ylim([0 60])
set(gca,'fontsize',16)
xlabel('Time (a.u.)')
ylabel('# Activated Receptors')
legend({sprintf('\\phi=%.2f',phi_vec(exPlot1)),sprintf('\\phi=%.2f',phi_vec(exPlot2))})
title('DiRT Simulations')
yticks(0:10:60)

subplot(1,2,2); hold on;
plot(tVec,alphaFn([JMax(exPlot1),tau_s(exPlot1)],tVec),...
    '--','LineWidth',1.5,'color','k')
plot(tVec,alphaFn([JMax(exPlot2),tau_s(exPlot2)],tVec),...
    '--','LineWidth',1.5,'color',colorScheme(1,:))
xlim([0 0.75]);
ylim([0 60])
title('Fitted Models')
set(gca,'fontsize',16)
xlabel('Time (a.u.)')
ylabel('# Activated Receptors')
legend({sprintf('\\phi=%.2f',phi_vec(exPlot1)),...
    sprintf('\\phi=%.2f',phi_vec(exPlot2))})
yticks(0:10:60)




%%
figure(4); clf; hold on;
plot(phi_vec,rsquared,'--','linewidth',1.5)
xlabel('\phi')
ylabel('R^2')
set(gca,'fontsize',16)
ylim([0 1])

