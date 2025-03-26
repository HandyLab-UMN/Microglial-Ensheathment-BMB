%%
% Make a raster plot 
% 1000 neurons are plotted, split between the population
%%
function plot_raster_fn(times, tinds, Ntot, Npop, Ncells, pinds, plot_num )

color_scheme = [38,30,101;
    0, 169, 69;
    255,172,16]/255;

color_scheme = repmat(color_scheme,ceil(Npop/3),1);
opacityLevel(1:3) = 1;
if Npop>3
    opacityLevel(4:6) = 0.5;
end

%%

s = [times; tinds];
s=s(:,s(1,:)>0); % eliminate wasted space in s

% Time window over which to plot in ms
% Tmin=890000; Tmax=901000;
Tmin=0; Tmax=5000;

% Plot spikes of 1000 neurons
if Ntot > 100000
    plot_per_total = 1000/Ntot; % percent of neurons to plot
else
    plot_per_total = 0.1;
end
index_start = zeros(Npop,1);
index_end = zeros(Npop,1);
n_per_pop = zeros(Npop,1);
for i = 1:Npop
    % Plot cells from E population
    n_per_pop(i) = Ncells(i)*plot_per_total;
%     n_per_pop(i) = 100;
    
    index_start(i) = pinds(i);
    index_end(i) = pinds(i)+n_per_pop(i);
    
    Iplot=find(s(1,:)>=Tmin & s(1,:)<=Tmax & s(2,:) > index_start(i) & s(2,:) <=index_end(i));
    
    % get the neuron indices
    neuroninds=s(2,Iplot);
    % map the neuron indices so that they lay next to each other
    if i > 1
        neuroninds = neuroninds - index_start(i) + sum(n_per_pop(1:(i-1)));
    end
    %neuroninds(neuroninds>Ne/2)=neuroninds(neuroninds>Ne/2)-Ne/2+nplot/2;
    
    if ~isempty(neuroninds)
        s_raster = double(s(1,Iplot))/1000;
        h(i) = scatter(s_raster,neuroninds,'Filled','SizeData',12,...
            'MarkerFaceColor', color_scheme(i,:),'MarkerEdgeColor', color_scheme(i,:),...
            'MarkerFaceAlpha',opacityLevel(i),'MarkerEdgeAlpha',opacityLevel(i));
        hold on;
    end
end

xlabel('Time (sec)')
yticklabels([])
set(gca,'fontsize',16)
    
if plot_num == 4
    if Npop == 4
        [lgd,icons] = legend([h(4), h(3), h(2), h(1)],{'VIP','SOM','PV','PN'},'fontsize',16);
    elseif Npop == 3
        if isempty(neuroninds)
            h(i) = plot(-1,-1,'.','MarkerSize',6,'color', color_scheme(i,:));
        end
        [lgd,icons] = legend([h(3), h(2), h(1)],{'SOM','PV','PN'},'fontsize',16);
    end

    if Npop>=3
        % Create a legend with 3 entries
        % Find the 'line' objects
        icons = findobj(icons,'Type','line');
        % % Find lines that use a marker
        icons = findobj(icons,'Marker','none','-xor');
        % % Resize the marker in the legend
        set(icons,'MarkerSize',20);
        % lgd.FontSize = 100;
    end
end

ylim([0 sum(n_per_pop)])

end



