%%
% Compares the firing rates and the power spectrum of the EIF simulations
% and the linear response theory
%
% Written by Gregory Handy, 09/18/2024
%%
clear; close all; clc;

%% Load parameter set
data_dir='./Data_sets/';

% Load the simulation data (trial 1)
% Trial 1 has awake, anesthetized, and emergence states
file_name = sprintf('EIF_Exp_Stats_Trial1.mat');
name_full = strcat(data_dir,file_name);
Sim_results1 = load(name_full);

% Load the simulation data (trial 2)
% Trial 2 has default, awake, anesthetized, and emergence states
% We ignore the default data for this plot
% Final result averages over these two trials
file_name = sprintf('EIF_Exp_Stats_Trial2.mat');
name_full = strcat(data_dir,file_name);
Sim_results2 = load(name_full);

% Load the linear response theory
file_name = sprintf('LR_Exp_Data.mat');
name_full = strcat(data_dir,file_name);
LR_results = load(name_full);

%%
num_conds = 4;
downsample = [1:3:80];

color_scheme = [38 30 101; 0 169 69; 255,172,16;...
    38 30 101; 0 169 69; 255,172,16]/255; 

f=figure(1); f.Position = [307 528 1077 387]; clf; 
clearvars h;
subplot(1,3,1); hold on;
for ii = 1:3

    rates_ave = (Sim_results1.rates_ave(ii,:) + Sim_results2.rates_ave(ii,2:end))/2;

    h(ii) = plot(1:3,rates_ave,...
        '.','markersize',16,'color',color_scheme(ii,:));
    plot(1:3,LR_results.rates(ii,2:end),...
        '--','markersize',16,'color',color_scheme(ii,:),'linewidth',1.5);
end
legend(h(1:3),{'E','PV','SST'},'location','northwest')
ylabel('Firing rate (Hz)')
xlabel('Brain state')
xticks([1 2 3])
set(gca,'fontsize',16)

percentChange = (LR_results.rates(:,end)-LR_results.rates(:,2))./LR_results.rates(:,2)*100;
fprintf('LR Percent change of exc. firing rates: %.2f\n',percentChange(1));
percentChange = (Sim_results1.rates_ave(:,end)-Sim_results1.rates_ave(:,1))./Sim_results1.rates_ave(:,1)*100;
fprintf('Sim Percent change of exc. firing rates: %.2f\n',percentChange(1));

% Plot the power spectrum
colorScheme = get(gca,'colororder');
colorScheme(3,:) = [0.5 0.5 0.5]; 
clearvars h legendString;
subplot(1,3,2); hold on;
for ii = 1:3
    Sim_ave_sfsp = (Sim_results1.sfsp1(ii,:) + Sim_results2.sfsp1(ii+1,:))/2;
    h(ii) = plot(Sim_results1.wfsp1(ii,downsample)/(2*pi)*1e3,Sim_ave_sfsp(1,downsample),...
        '.','markersize',16,'color',colorScheme(ii,:));

    plot(LR_results.params.omega*1e3,real(squeeze(LR_results.yy_freq(1,1,:,ii+1))),...
        '--','linewidth',1.5,'color',colorScheme(ii,:));
    xlim([0 100])
    set(gca,'fontsize',16)
end
legend(h,{'Awake','Anesthetized','Emergence'})
xlabel('Freq (Hz)')
ylabel('C_{ee}(f) (spikes^2 per Hz)')


% Plot the coherence
clearvars h legendString;
subplot(1,3,3); hold on;
for ii = 1:3

    % Sim coherence
    powerSpecE1_Sim1 = Sim_results1.sfsp1(ii,:);
    powerSpecE2_Sim1 = Sim_results1.sfsp2(ii,:);
    crossSpec_Sim1 = Sim_results1.sfsp12(ii,:);
    powerSpecE1_Sim2 = Sim_results2.sfsp1(ii+1,:);
    powerSpecE2_Sim2 = Sim_results2.sfsp2(ii+1,:);
    crossSpec_Sim2 = Sim_results2.sfsp12(ii+1,:);

    
    powerSpecE1_ave = (powerSpecE1_Sim1 + powerSpecE1_Sim2)/2;
    powerSpecE2_ave = (powerSpecE2_Sim1 + powerSpecE2_Sim2)/2;
    crossSpec_ave = (crossSpec_Sim1 + crossSpec_Sim2)/2;
    

    coher_Sim = abs(crossSpec_ave).^2./(powerSpecE1_ave.*powerSpecE2_ave);

    % LR coherence
    powerSpecE1_LR = squeeze(LR_results.yy_freq(1,1,:,ii+1));
    powerSpecE2_LR = squeeze(LR_results.yy_freq(4,4,:,ii+1));
    crossSpecE2_LR = squeeze(LR_results.yy_freq(1,4,:,ii+1));
    coher_LR = real(abs(crossSpecE2_LR).^2./(powerSpecE1_LR.*powerSpecE2_LR));

    h(ii) = plot(Sim_results1.wfsp1(ii,downsample)/(2*pi)*1e3,coher_Sim(1,downsample),...
        '.','markersize',16,'color',colorScheme(ii,:));
    plot(LR_results.params.omega*1e3,coher_LR,...
        '--','linewidth',1.5,'color',colorScheme(ii,:));
    xlim([0 100])
    set(gca,'fontsize',16)
end
legend(h,{'Awake','Anesthetized','Emergence'})
xlabel('Freq (Hz)')
ylabel('c_{ee}(f) (spikes^2 per Hz)')







