%%
% Finds the cross-correlation function (corrAO_sim) and power spectrum 
% (sfsp) of the neuron population that lies between the indices 
% index_state and index_end
%
% Creates a "master" time series for the neuron population contains the
% spike counts in 1 msec time bins
%
% If find_adj == 1, the auto-corr is subtracted from the xcorr of this 
% master time seris  
%   Note: A big-time sink is finding the auto-correlations for each 
%   spike train
%
% Can be made faster, but less accurate by setting num_sampled < total_neurons
%
% Note: wfsp (freq.) has units of radians/msec
%%
function [corrAO_adj_sim, corrAO_sim, lags, sfsp_adj, sfsp, wfsp] = ...
    power_spec_fn(spkIndex,spkTime, index_start, index_end, total_neurons, ...
    num_sampled,T_start,T_end,find_adj)

%% shrink the spike index and time by keeping the nonzero ones
spkIndex = spkIndex(spkTime>0)+1;
spkTime = spkTime(spkTime>0);

spkIndex_corr = spkIndex;
spkTime_corr = spkTime;

%%
bin_corr = 1; %ms; changing this is not recommended (need to rescale other numbers throughout)
tEdge = T_start : bin_corr : T_end;
neuronEdge = [index_start-0.5:1:index_end+0.5];

%% binned spike count
bSpk = histcounts2(spkIndex_corr, spkTime_corr, neuronEdge, tEdge);

%% Restrict the estimate to the number of sampled cells
if num_sampled < total_neurons
    % Find those neurons with a non-zero firing rate
    Igood=find(rates(index_start:index_end)>=0.1);
    
    temp=randperm(numel(Igood),num_sampled);
    Inds=Igood(temp);
    
    % Keep the "good" spike counts
    bSpk=bSpk(Inds,:);
end

%%
tMax = 250;
tMax_bin = tMax/bin_corr;
% minus the mean for spk train of population 1
sum_spkTrain = sum(bSpk,1);
sum_spkTrain = sum_spkTrain - mean(sum_spkTrain);

[total_corr, lags] = xcorr(sum_spkTrain, tMax_bin); 

%% Find the auto-correlation of each spike train
if find_adj == 1
    auto_corr = 0;
    parfor kk = 1:num_sampled
        single_spkTrain = bSpk(kk,:)-mean(bSpk(kk,:));
        [temp_corr,~] = xcorr(single_spkTrain, single_spkTrain, tMax_bin);
        auto_corr = auto_corr + temp_corr;
    end
    
    %% Subtract the auto-corr from the total-corr
    cross_corr_adj = total_corr - auto_corr;
    cross_corr_adj_normalized = cross_corr_adj / ((T_end-T_start)*(num_sampled)*(num_sampled-1));
    corrAO_adj_sim = cross_corr_adj_normalized;
else
    corrAO_adj_sim = nan; % set just so the function has something to return
end

%% 
corrAO_sim = total_corr / ((T_end-T_start)*num_sampled^2);


%% Take the Fourier transform
fmaxHz=100;
T_lim=400;

ft_index_start = find(lags==-100,1);
ft_index_end = find(lags==100,1);

nffsp=floor(2*T_lim*fmaxHz/1000); 
ffsp=(1/(2*T_lim))*[1:nffsp]';
wfsp=ffsp*2*pi;
sfsp_adj=zeros(nffsp,1);
sfsp=zeros(nffsp,1);
for k=1:nffsp
    if find_adj == 1
        sfsp_adj(k) = trapz(lags(ft_index_start:ft_index_end),...
            exp(-1i*wfsp(k)*lags(ft_index_start:ft_index_end)).*...
            corrAO_adj_sim(ft_index_start:ft_index_end));
    end

    sfsp(k) = trapz(lags(ft_index_start:ft_index_end),...
        exp(-1i*wfsp(k)*lags(ft_index_start:ft_index_end)).*...
        corrAO_sim(ft_index_start:ft_index_end));
end

sfsp_adj = real(sfsp_adj); % cleans up and imaginary bits
sfsp = real(sfsp); % cleans up and imaginary bits

end

