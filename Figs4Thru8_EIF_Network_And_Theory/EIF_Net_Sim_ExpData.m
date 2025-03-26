%%
% This runs the spiking network simulation using mex code.
%
% Saves the full spike train of all trials, and then condenses the data to
% find the average steady state firing rate, the power-spectrum and the
% cross-spectrum
%
% Denote the desired location to save the large data file (i.e., all spikes)
% in data_dir
%
% Written by Gregory Handy, 07/26/2024
%%
clear; clc; %close all;

restoredefaultpath;
folder = fileparts(which('EIF_Net_Sim_ExpData.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Compile the mex code

% Use to compile the mex code (do after each edit to the mex file)
% Run the following if running mex for the first time
% run mex -setup -c
mex ./Mex_Functions/EIF_Mex_AlphaFnPlusEns.c

%% Save the data in this directory
% Warning: Can be a large file!
data_dir='./Data_sets/';

%% Load the parameters and seed the random number generator
surround_diff = 0; % 0 (iso), 90 (cross)
ffwdInputStr = [0.75 0.75 0 0.75]; % Strength of the feedforward inputs
r_vip = 6; % VIP default firing rate
params = EIF_params_official_fn(surround_diff,0.75,{'range',r_vip},'med');

% Seed the random number generator
rng('shuffle');
s = rng();
random_seed = s.Seed;

%% Ensheathment parameters

% Columns are the probability of being/no being ensheathed
%   1st column: probability of not being ensheathed
%   2nd column: probability of being ensheathed
% Rows are the types of synapses being ensheathed
%   1st: Exc synapses
%   2nd row: PV synapses
%   3rd row: SST synapses
% Note: all rows should add up to 1
prob_en_tensor = zeros(3,4,3);

% Default/no ensheathment
prob_en_tensor(:,:,1) = [1.00 0.00 0.00 0.00;
                         1.00 0.00 0.00 0.00;
                         1.00 0.00 0.00 0.00];

% Awake
normFactor = sum([0.42 0.14 0.06]);
prob_en_tensor(:,:,2) = [1.00 0.00 0.00 0.00;
                         0.80 0.2*0.42/normFactor .2*0.14/normFactor .2*0.06/normFactor;
                         0.80 0.2*0.42/normFactor .2*0.14/normFactor .2*0.06/normFactor;];

% Anesthesia
normFactor = sum([0.42 0.30 0.16]);
prob_en_tensor(:,:,3) = [1.00 0.00 0.00 0.00;
                         0.20 0.8*0.42/normFactor 0.8*0.30/normFactor 0.8*0.16/normFactor;
                         0.20 0.8*0.42/normFactor 0.8*0.30/normFactor 0.8*0.16/normFactor];

% Emergence
normFactor = sum([0.49 0.23 0.11]);
prob_en_tensor(:,:,4) = [1.00 0.00 0.00 0.00;
                         0.267 0.49*.733/normFactor 0.23*.733/normFactor 0.11*.733/normFactor;
                         0.267 0.49*.733/normFactor 0.23*.733/normFactor 0.11*.733/normFactor];

% Repeat these parameter values for the surround populations
prob_en_tensor = repmat(prob_en_tensor,params.Npop/3,1);
if sum(sum(abs(sum(prob_en_tensor,2)-1)>1e-10)) ~=0
    error('Probabilities do not all add up to 1');
end

% num_conds = length(ffwdInputStr);
num_conds = 4;

%% Run the simulations
tic;
fprintf('Running the spiking simulations \n');
parfor cond_num = 1:num_conds

    params = EIF_params_official_fn(surround_diff,ffwdInputStr(cond_num),{'range',r_vip},'med');

    % Strength of ensheathment
    params.s_en_vec = [0.0 0.33 0.67 1.0]; % Can be adjusted
    s_en = params.s_en_vec;
    prob_en = prob_en_tensor(:,:,cond_num);

    %% Create the connectivity matrix to be used for all simulations

    fprintf('Creating the connectivity matrix \n');
    [wind,wipost,wstr,ensh_realizations,tausyn] = ...
        Diff_approx_gen_weights_ens(params.Ncells,params.p,params.J,...
        prob_en,s_en,params.tau_s,params.dt);

    % offset the index by one for the MEX code (C starts indexing at 0)
    wind_mex = int32(wind-1);
    wipost_mex = int32(wipost-1);
    pinds_mex = int32(params.pinds-1);

    %% Define the feedforward input (mean and variance)
    mu_ext = params.mu_bg + params.mu_stim + params.mu_vip;
    sigma_ext = sqrt(params.sigma_bg.^2+params.sigma_stim.^2+...
        params.sigma_vip.^2).*sqrt(2*params.tau_m);

    %% Run the simulation
    [rates_trial,times,tinds] = ...
        EIF_Mex_AlphaFnPlusEns(params.T, params.NT, params.recstart,...
        params.Npop, params.Ntot, params.Ncells, params.dt, params.rates, params.p, params.J,...
        tausyn, params.tau_m, params.EL, params.vth, params.vre, params.tau_ref,...
        wind_mex, wipost_mex, wstr, pinds_mex,random_seed,mu_ext,...
        sigma_ext,params.DeltaT,params.VT, params.tauDelay,...
        params.maxDelay,params.shared_noise_mag.*sqrt(2*params.tau_m),...
        params.maxns,length(s_en),ensh_realizations);

    % Save the spike times and spike indices
    spkMatrix = [times(times>0);tinds(times>0)];

    save_data_exp(data_dir,random_seed,rates_trial,spkMatrix,params,cond_num)
end
toc;

%% Find the steady state firing rates and plot the results

fprintf('Loading the data and estimating firing rates \n')
rates_ave = zeros(params.Npop,num_conds);
for ii = 1:num_conds
    % Load the data
    file_name = sprintf('EIF_ExpData_%d_%d',ii,random_seed);
    name_full = strcat(data_dir,file_name);
    load(name_full,'rates_trial')
    % Estimate the firing rate
    for jj = 1:params.Npop
        rates_ave(jj,ii) = mean(rates_trial(1,params.pinds(jj):params.pinds(jj+1)-1));
    end
end

color_scheme = [38 30 101; 0 169 69; 255,172,16]/255;
color_scheme = repmat(color_scheme,params.Npop/3,1);
clearvars h;
figure(1); clf; hold on;
for ii = 1:3
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

set(gca,'fontsize',16)
legend(h(1:3),{'E','PV','SOM'})
ylabel('Firing rate (Hz)')
% xlabel('s_{en}')
% xlim([min(params.s_en_vec) max(params.s_en_vec)])
% title(sprintf('Prob PV Ensheathment: %.2f',params.prob_en(2,2)))

%% Plot an example raster plot from this simulation

figure(77); clf;
file_name = sprintf('EIF_ExpData_%d_%d',1,random_seed);
name_full = strcat(data_dir,file_name);
load(name_full,'rates_trial','spkMatrix','params')
plot_raster_fn(spkMatrix(1,:), spkMatrix(2,:), params.Ntot, 3,...
    params.Ncells, params.pinds, 1)
xlim([4 4.5])


