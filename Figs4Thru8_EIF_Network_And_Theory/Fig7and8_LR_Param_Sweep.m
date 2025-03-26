%%
% Plot the results from the parameter sweeps for the two data sets
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('Fig7and8_LR_Param_Sweep.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Load each of the data sets and create the plots
ensPop = {'PV','SST'};
for kk = 1:2

    %% Load the data sets
    if kk == 1
        legendText = {'\rho_p'};
    else
        legendText = {'\rho_s'};
    end

    data_dir='./Data_sets/';
    file_name = sprintf('LR_%s_Data.mat',ensPop{kk});
    name_full = strcat(data_dir,file_name);
    load(name_full)

    %% Plot the population firing rates
    colorScheme = get(gca,'ColorOrder');

    num_p_conds = length(prob_en_vec);
    f1 = figure(1+2*(kk-1));
    f1.Position =  [304 368 392 596];
    for ii = 1:num_p_conds
        for jj = 1:3
            subplot(3,1,mod(jj+2,3)+1); hold on;
            h(ii,jj) = plot(params.s_en_vec,rates(jj,:,ii),...
                '-','color',[colorScheme(mod(jj+2,3)+1,:) ii/num_p_conds],...
                'linewidth',1.5);
            set(gca,'fontsize',16)
            xlabel('s_{en}')
        end
    end

    for jj = 1:3
        subplot(3,1,jj);

        if jj == 1
            ylabel('Exc. Firing rate (Hz)')
            legend(h([1 end],jj),{strcat(legendText{1},'=0'),strcat(legendText{1},'=0.8')});
            ylim([4 10])
        elseif jj == 2
            ylabel('PV Firing rate (Hz)')
            ylim([5 10])
        else
            ylabel('SST Firing rate (Hz)')
            ylim([6 10])
        end
    end


    %%  Estimate and plot the gamma power
    pos_indices = find(params.omega*1e3>10);
    [gamma_power, gamma_index] = ...
        max(real(squeeze(yy_freq(1,1,pos_indices,:,:))));
    gamma_power = squeeze(gamma_power);
    gamma_index = squeeze(gamma_index);

    f2 = figure(2+2*(kk-1)); f2.Position = [697 373 713 591];
    % Plot the gamma power examples
    subplot(2,2,1); hold on;
    plot(params.omega*1e3,real(squeeze(yy_freq(1,1,:,1,1))),'linewidth',1.5,'color','k')
    plot(params.omega*1e3,real(squeeze(yy_freq(1,1,:,8,6))),'--','linewidth',1.5,'color','k')
    xlim([0 100])
    % Plotting the points corresponding to (params.s_en_vec(8), prob_en_vec(6))
    legend(strcat(legendText{1},'=0'),strcat('s_{en} = 0.7, ',legendText{1},'=0.5'))
    set(gca,'fontsize',16)
    ylabel('c_{ee}(f) (spikes^2 per Hz)')
    xlabel('Freq (Hz)')
    legend AutoUpdate off
    plot(params.omega(pos_indices(gamma_index(8,6)))*1e3,gamma_power(8,6),'k*',...
        'markersize',6,'LineWidth',1.5)
    ylim([0 2.5e-5])

    % Plot the gamma power over prob x s_en space
    [prob,s_en] = meshgrid(prob_en_vec,params.s_en_vec);
    subplot(2,2,2); hold on;
    surf(prob,s_en,(gamma_power-gamma_power(1,1))./gamma_power(1,1)*100,'FaceColor','interp')
    colorbar
    clim([0 400])
    colormap(linspecer)
    xlim([0.2 0.8])
    ylim([0 1])
    xlabel(legendText{1})
    ylabel('s_{en}')
    set(gca,'fontsize',16)
    title('Gamma power')
    plot3(params.s_en_vec(6),prob_en_vec(8),...
        max(max((gamma_power-gamma_power(1,1))./gamma_power(1,1)*100))*2,...
        'k*','markersize',8,'LineWidth',1.5)

    %% Calculate and plot the coherence
    powerSpecE1 = real(squeeze(yy_freq(1,1,:,:,:)));
    powerSpecE2 = real(squeeze(yy_freq(4,4,:,:,:)));
    crossSpec = squeeze(yy_freq(1,4,:,:,:));
    coher = abs(crossSpec).^2./(powerSpecE1.*powerSpecE2);
    max_coher = zeros(size(gamma_index));
    for ii = 1:num_p_conds
        for jj = 1:length(params.s_en_vec)
            max_coher(jj,ii) = coher(pos_indices(gamma_index(jj,ii)),jj,ii);
        end
    end

    % Plot the coherence examples
    subplot(2,2,3); hold on;
    plot(params.omega*1e3,real(squeeze(coher(:,1,1))),'linewidth',1.5,'color','k')
    plot(params.omega*1e3,real(squeeze(coher(:,8,6))),'--','linewidth',1.5,'color','k')
    plot(params.omega(pos_indices(gamma_index(8,6)))*1e3,...
        real(squeeze(coher(pos_indices(gamma_index(8,6)),8,6))),'.','markersize',16,'color','k')
    xlim([0 100])
    set(gca,'fontsize',16)
    ylabel('Coherence')
    xlabel('Freq (Hz)')
    ylim([0 1])

    % Plot the coherence over prob x s_en space
    ax = subplot(2,2,4); hold on;
    cMin = 0; cMax = 50;
    gColorScheme = gBlendFn(cMin,cMax);
    max_coher_default = max_coher(1,1);
    surf(prob,s_en,(max_coher-max_coher_default)/max_coher_default*100,'FaceColor','interp')
    colorbar
    clim([cMin cMax])
    colormap(ax, gColorScheme)

    plot3(params.s_en_vec(6),prob_en_vec(8),max(max((max_coher-max_coher_default)/max_coher_default*100)),...
        'k.','markersize',20,'LineWidth',1.5)

    xlim([0.2 0.8])
    ylim([0 1])
    xlabel(legendText{1})
    ylabel('s_{en}')
    title('Gamma coherence')
    set(gca,'fontsize',16)
end

%% To save this panels, we need to render it with 'painters'
% set(gcf, 'renderer','painters');
% figure(1); print('-vector','PV_rates',"-dpdf")
% figure(2); print('-vector','PV_corr',"-dpdf")
% figure(3); print('-vector','SST_rates',"-dpdf")
% figure(4); print('-vector','SST_corr',"-dpdf")


