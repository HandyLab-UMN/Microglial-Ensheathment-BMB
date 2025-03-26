%%
% Runs the linear response mean-field approximation to find the
% population average firing rates, power spectrum and coherence
% using ensheathment parameters based on the experimental data 
%
% Written by Gregory Handy, 09/18/2024
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('Linear_Resp_Calc_ExpData.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Compile mex code (if necessary)

mex ./Mex_functions/calc_Rate_cuEIF.cpp
mex ./Mex_functions/calc_Power_cuEIF.cpp
mex ./Mex_functions/calc_Susc_cuEIF.cpp

%% Load parameters
surround_diff = 0; % 0 (iso), 90 (cross)
ffwdInputStr = [0.75 0.75 0 0.75]; % Strength of the feedforward inputs
r_vip = 6; % VIP default firing rate
params = EIF_params_official_fn(surround_diff,0,{'range',r_vip},'med');

%% Ensheathment parameters
% Columns are the probability of being/no being ensheathed
%   1st column: probability of not being ensheathed
%   2nd column: probability of being ensheathed
% Rows are the types of synapses being ensheathed
%   1st: Exc synapses
%   2nd row: PV synapses
%   3rd row: SST synapses
% Note: all rows should add up to 1
prob_en_tensor = zeros(3,4,4);

% Default/No ensheathment
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
if sum(sum(sum(prob_en_tensor,2)~=1)) ~=0
    error('Probabilities do not all add up to 1');
end

num_conds = length(ffwdInputStr);

%% Preallocate 
rates = zeros(params.Npop,num_conds);
yy_freq = zeros(params.Npop,params.Npop,params.bins,num_conds);
yy_time = zeros(params.Npop,params.Npop,params.bins,num_conds);

%% Loop over the stimuli
parfor jj = 1:num_conds
    
    params = EIF_params_official_fn(surround_diff,ffwdInputStr(jj),{'range',r_vip},'med');

    % Strength of ensheathment
    params.s_en = [0.0 0.33 0.67 1.0];
    
    % Probability of ensheathment
    params.prob_en = prob_en_tensor(:,:,jj);
    
    % Compute the linear response theory
    [rates(:,jj),yy_freq(:,:,:,jj),yy_time(:,:,:,jj)]=...
        LinearResponse_fn(params,params.s_en,params.prob_en,1);
end

%% Save the results
data_dir='./Data_sets/';
file_name = sprintf('LR_Exp_Data_NEW.mat');
name_full = strcat(data_dir,file_name);
save(name_full,'rates','yy_freq','yy_time','params')

%% Plot the steady state firing rates
color_scheme = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880;...
    0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880]; 

f = figure(1); clf; f.Position = [532 386 1110 417];
subplot(1,3,1); hold on;
for ii = 1:params.Npop
    h(ii) = plot(1:num_conds, rates(ii,:),'--s','markersize',8,...
        'linewidth',1.5,'color',color_scheme(ii,:));
end
set(gca,'fontsize',16)
legend(h(1:3),{'E','PV','SOM'})
xticks(1:4)
xticklabels({'No Ensh.','Awake','Anesthetized','Emergence'})
ylabel('Firing rate (Hz)')

colorScheme = [0.4940    0.1840    0.5560; 0 0.447 0.7410; 0.85 0.325 0.098; 0.5 0.5 0.5];
% Plot the power spectrums
subplot(1,3,2); hold on;
for kk = 1:num_conds
    plot(params.omega*1e3,(real(squeeze(yy_freq(1,1,:,kk)))),...
        '-','linewidth',1.5,'color',colorScheme(kk,:))
    xlim([0 100])

    pos_indices = find(params.omega*1e3>0);
    [val, index] = max(real(squeeze(yy_freq(1,1,pos_indices,kk))));
    max_freq(kk) = params.omega(pos_indices(index))*1e3;
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{ee}(f) (spikes^2 per Hz)')
legend({'No Ensh.','Awake','Anesthetized','Emergence'})

subplot(1,3,3); hold on;
clearvars h;
for kk = 1:num_conds
    powerSpecE1_LR = squeeze(yy_freq(1,1,:,kk));
    powerSpecE2_LR = squeeze(yy_freq(4,4,:,kk));
    crossSpecE2_LR = squeeze(yy_freq(1,4,:,kk));
    coher_LR = real(abs(crossSpecE2_LR).^2./(powerSpecE1_LR.*powerSpecE2_LR));

    h(kk) = plot(params.omega*1e3,coher_LR,...
        '-','linewidth',1.5,'color',colorScheme(kk,:));
    xlim([0 100])
    set(gca,'fontsize',16)
end

ylabel('coherence')
xlabel('Freq (Hz)')
set(gca,'fontsize',16)
legend({'No Ensh.','Awake','Anesthetized','Emergence'})
ylim([0 1])


fprintf('Percent change in exc. firing ratesL %.2f%%\n',(rates(1,4)-rates(1,2))./rates(1,2)*100)

