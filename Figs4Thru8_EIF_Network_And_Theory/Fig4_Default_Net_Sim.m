%%
% Plot the results from the default network (i.e., no ensheathment)
% Loads the data and computers the averaged firing rates over rolling
% windows
%
% Written by Gregory Handy, 09/18/2024
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('Fig4_Default_Net_Sim.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Load the data

% Data directory path
data_dir = './Data_sets/';

% Load the simulation data
file_name = sprintf('EIF_Exp_Spikes_Default.mat');
name_full = strcat(data_dir,file_name);
Sim_results = load(name_full);

% Load the corresponding statistics from this simulation
% This trial has the default case of no ensheathment
file_name = 'EIF_Exp_Stats_Trial2'; 
name_full = strcat(data_dir,file_name);
Sim_stats = load(name_full);

%%
dataSet = Sim_results;
totalExamples = 1;
sliding_window = 10;
tvec = [3900:sliding_window:5000];
reSim_roll1 = zeros(length(tvec),1);
reSim_roll2 = zeros(length(tvec),1);
for tt = 1:length(tvec)
    reSim_roll1(tt) = 1000*nnz(dataSet.spkMatrix(1,:)>sliding_window*(tt-1) & ...
        dataSet.spkMatrix(1,:)<sliding_window*tt & ...
        dataSet.spkMatrix(2,:)<=dataSet.params.Ncells(1))/(dataSet.params.Ncells(1)*(sliding_window));
    
    reSim_roll2(tt) = 1000*nnz(dataSet.spkMatrix(1,:)>sliding_window*(tt-1) & ...
        dataSet.spkMatrix(1,:)<sliding_window*tt & ...
        dataSet.spkMatrix(2,:)>=5000 & dataSet.spkMatrix(2,:)<9000)/(dataSet.params.Ncells(1)*(sliding_window));
end


%%
f = figure(1); f.Position = [870 241 489 675]; clf;

subplot(3,1,1); hold on;
plot(tvec/1e3, reSim_roll1,'linewidth',1.5,'color',[38,30,101]/255)
plot(tvec/1e3, reSim_roll2,'linewidth',1.5,'color',[38/255,30/255,101/255,0.5])
xlim([4 4.5])
xticks(4:.1:4.5)
xticklabels(0:.1:1)
ylim([2 9])
set(gca,'fontsize',16)
ylabel('Firing Rate (Hz)')
legend('Exc. Pop 1','Exc. Pop 2')
box off

subplot(3,1,2:3); hold on;
plot_raster_fn(dataSet.spkMatrix(1,:),dataSet.spkMatrix(2,:),...
        dataSet.params.Ntot,6,dataSet.params.Ncells,...
        dataSet.params.pinds,1)
xlim([4 4.5])
ylabel('Neuron Index')
xticks(4:.1:4.5)
xticklabels(0:.1:1)
% ylim([500 1000])
box off

%%
colorScheme = get(gca,'colororder');
f = figure(2); f.Position = [544 530 1112 420]; clf; 

% Auto-correlation function
subplot(1,3,1); hold on;
plot(Sim_stats.lags1(1,:), Sim_stats.corrAO_sim1(1,:),'-','linewidth',1.5,'color',colorScheme(1,:))
set(gca,'fontsize',16)
xlim([-100 100])
xticks([-100 0 100])
xticklabels([0.1 0 0.1])
ylim([-0.05*1e-5 0.3*1e-5])
xlabel('Time lag (seconds)')
ylabel('C_{e_1e_1}(t) (Hz^2)')

% Sim power and cross-spectral
subplot(1,3,2); hold on;
plot(Sim_stats.wfsp1(1,:)/(2*pi)*1e3, Sim_stats.sfsp1(1,:),'-','linewidth',1.5,'color',colorScheme(1,:))
plot(Sim_stats.wfsp1(1,:)/(2*pi)*1e3, real(Sim_stats.sfsp12(1,:)),'k-','linewidth',1.5)
legend('Power Spectrum','Cross Spectrum')
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('C_{e_ie_j}(f) (spikes^2/second)')
ylim([0 0.15*1e-4])


subplot(1,3,3); hold on;
% Calculate sim coherence
powerSpecE1_Sim = Sim_stats.sfsp1(1,:);
powerSpecE2_Sim = Sim_stats.sfsp2(1,:);
crossSpec_Sim = Sim_stats.sfsp12(1,:);
coher_Sim = abs(crossSpec_Sim).^2./(powerSpecE1_Sim.*powerSpecE2_Sim);
plot(Sim_stats.wfsp1(1,:)/(2*pi)*1e3,coher_Sim,...
        '-','linewidth',1.5,'color',colorScheme(1,:));
ylim([0 1])
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('Coherence')




