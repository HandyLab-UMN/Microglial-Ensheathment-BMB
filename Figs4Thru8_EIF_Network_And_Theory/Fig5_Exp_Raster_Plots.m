clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('Fig5_Exp_Raster_Plots.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Load the three data sets to creat the raster plots

expStates = {'Awake','Anesthetized','Emergence'};

fprintf('Loading the three data sets\n');
data_dir = './Data_sets/';
for ii = 1:length(expStates)
    file_name = sprintf('EIF_Exp_Spikes_%s',expStates{ii});
    name_full = strcat(data_dir,file_name);
    expSims{ii} = load(name_full);
end
%%
tic;
f = figure(1); f.Position = [400 223 1233 749]; clf; hold on;

for ii = 1:length(expStates)

    dataSet = expSims{ii};

    subplot(2,4,ii+5)
    plot_raster_fn(dataSet.spkMatrix(1,:),dataSet.spkMatrix(2,:),...
        dataSet.params.Ntot,1,dataSet.params.Ncells,...
        dataSet.params.pinds,1)
    ylabel('Neuron Index')
    xlim([4 4.25])
    xlim([4 4.25])
    xticks([4 4.1 4.2])
    xticklabels({0,0.1,0.2})
    box off
end
toc;

%% Get the rolling averages for the firing rates
sliding_window = 25; % length of the averaging window
tvec = 3900:sliding_window:4500; % time frame for the averaging
reSim_roll = zeros(length(tvec),length(expStates));
for ii = 1:length(expStates)
    dataSet = expSims{ii};

    for tt = 1:length(tvec)
        reSim_roll(tt,ii) = 1000*nnz(dataSet.spkMatrix(1,:)>sliding_window*(tt-1) & ...
            dataSet.spkMatrix(1,:)<sliding_window*tt & ...
            dataSet.spkMatrix(2,:)<=dataSet.params.Ncells(1))/(dataSet.params.Ncells(1)*(sliding_window));
    end

    subplot(2,4,ii+1)
    plot(tvec/1000,reSim_roll(:,ii),'linewidth',1.5,'color',[38,30,101]/255)
    ylabel('Firing Rates (Hz)')
    xlabel('Time (sec)')
    xlim([4 4.25])
    xticks([4 4.1 4.2])
    xticklabels({0,0.1,0.2})
    ylim([-0.5 13])
    set(gca,'fontsize',16)
    box off
    title(expStates{ii})
end

%% Recreate the modified experimental figures
subplot(2,4,1)
bar([20 80; 80 20; 73.3 26.7],'stacked')
set(gca,'fontsize',16)
yticks([0:20:100])
box off
legend('High contact','No/low contact')
ylabel('Neuron (%)')

subplot(2,4,5)

bar([38 12, 17; 42 42 48; 14 30 23; 6 16 11]);
set(gca,'fontsize',16)
legend('Awake', 'Isoflurane', 'Emergence')
xticklabels({'No contact','Contact','Enwrapped','Engulfed'})
box off

