%%
% Estimates the firing rates and power spectrum of the exc. neurons
% Saves the power spectrum data
%
% Note: For the power spectrum to be at least a little accurate, 
%       params.T needs to be at least 1e6
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('Analyze_Spike_Data.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Load parameter values

data_dir = 'Data_sets/';
% data_dir = '/Users/ghandy/Local_Data/Raw_Spike_Trains/';
% random_seed = 609012141;
random_seed = 1057695155;
% random_seed = 607333488;

file_name = sprintf('EIF_ExpData_%d_%d',1,random_seed);
name_full = strcat(data_dir,file_name);
load(name_full,'params')

%% Find the steady state firing rates
num_conds = 4;
rates_ave = zeros(params.Npop,num_conds);
cov_est = zeros(num_conds,1);
% Estimate the firing rate
for ii = 1:num_conds

    fprintf('Estimating firing rates and long-time cov for dataset %d\n',ii);

    % Load the data    
    file_name = sprintf('EIF_ExpData_%d_%d',ii,random_seed);
    name_full = strcat(data_dir,file_name);
    load(name_full,'rates_trial','spkMatrix')
    
    for jj = 1:params.Npop
        rates_ave(jj,ii) = mean(rates_trial(1,params.pinds(jj):params.pinds(jj+1)-1));
    end

    % Find the long-time covariance (for comparison at freq = 0)
    [~, scc_all, ~, ~]=...
        get_statistics(spkMatrix,params.T,params.pinds,...
        params.Ncells,params.recstart);
    cov_est(ii) = mean(scc_all)/250;
end

%% Find the power-spectrum 
% Note: Requires a decent amount of memory
% May take a little to run
clearvars corrAO_sim1 lags1 sfsp1 wfsp1
for mod_num = 1:num_conds
    fprintf('Finding the power-spectrum of the E populations for dataset %d\n',mod_num);
    file_name = sprintf('EIF_ExpData_%d_%d',mod_num,random_seed);
    name_full = strcat(data_dir,file_name);
    sim_mod = load(name_full,'rates_trial','spkMatrix','params');
    sim_mod.times = sim_mod.spkMatrix(1,:);
    sim_mod.tinds = sim_mod.spkMatrix(2,:);

    [~, corrAO_sim1(mod_num,:), lags1(mod_num,:), ~, sfsp1(mod_num,:), wfsp1(mod_num,:)] = ...
        power_spec_fn(sim_mod.tinds,sim_mod.times, sim_mod.params.pinds(1),...
        sim_mod.params.pinds(2)-1, sim_mod.params.Ncells(1), sim_mod.params.Ncells(1), ...
        params.recstart,params.T,0);


    fprintf('Estimating power-spectrum for for the surround E Pop: Stim %d\n',mod_num);
    [~, corrAO_sim2(mod_num,:), lags2(mod_num,:), ~, sfsp2(mod_num,:), wfsp2(mod_num,:)] = ...
        power_spec_fn(sim_mod.tinds,sim_mod.times, sim_mod.params.pinds(4),...
        sim_mod.params.pinds(5)-1, sim_mod.params.Ncells(4), sim_mod.params.Ncells(4), ...
        sim_mod.params.recstart,sim_mod.params.T,0);


    fprintf('Estimating cross-spectrum: Stim %d \n',mod_num);
    [corrAO_sim12(mod_num,:),lags12(mod_num,:), sfsp12(mod_num,:), wfsp12(mod_num,:)] = cross_spectrum_fn(sim_mod.tinds,sim_mod.times,...
        sim_mod.params.pinds(1), sim_mod.params.pinds(2)-1, sim_mod.params.pinds(4),...
        sim_mod.params.pinds(5)-1,sim_mod.params.recstart,sim_mod.params.T,...
        sim_mod.params.Ncells(1),sim_mod.params.Ncells(4));
end

file_name = sprintf('./Data_sets/EIF_Exp_Data_%d.mat',random_seed);
% save(file_name,'corrAO_sim1','lags1','wfsp1','sfsp1',...
%     'rates_ave','params','cov_est')
save(file_name,'corrAO_sim1','lags1','wfsp1','sfsp1',...
    'corrAO_sim2','lags2','wfsp2','sfsp2',...
    'corrAO_sim12','lags12','wfsp12','sfsp12',...
    'rates_ave','params','cov_est')


%% Plot an example power spectrum for s_en
figure(1); %clf;
subplot(1,2,1); hold on;
plot(lags1(1,:),corrAO_sim1(1,:),'linewidth',1.5)
xlabel('Lags (msec)')
set(gca,'fontsize',16)
ylabel('Autocorrelogram')

subplot(1,2,2); hold on;
plot(wfsp1(1,:)/(2*pi)*1e3,sfsp1(4,:),'.','markersize',16);
plot(0,cov_est(1),'*','markersize',16,'linewidth',1.5)

legend('Simulation','Long time covariance est')
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{e1e1}(f)')

%% Plot the firing rates
color_scheme = [38 30 101; 0 169 69; 255,172,16]/255;
color_scheme = repmat(color_scheme,params.Npop/3,1);
clearvars h;
figure(11); hold on;
for ii = 1:params.Npop
    if ii <= 3
        markerType = '.';
        markerSize = 16;
    else
        markerType = 's';
        markerSize = 8;
    end
    h(ii) = plot(1:num_conds,rates_ave(ii,:),markerType,'MarkerFaceColor',color_scheme(ii,:),...
        'markersize',markerSize,'color',color_scheme(ii,:));
end

legend(h(1:3),{'E','PV','SST'})
ylabel('Firing rate (Hz)')
xlabel('s_{en}')
set(gca,'fontsize',16)

