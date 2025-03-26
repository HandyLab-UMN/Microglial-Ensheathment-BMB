%%
% Finds the cross-correlation (corrAO_sim) and cross-spectrum (sfsp)
% across the two neuron populations
%
% Creates a "master" time series for each population that contains the
% spike counts in 1 msec time bins for all cells in each population
%
% Note: wfsp (freq.) has units of radians/msec
%%
function [corrAO_sim,lags, sfsp, wfsp] = cross_spectrum_fn(spkIndex,spkTime,...
    index_start_pop1, index_end_pop1, index_start_pop2, index_end_pop2,...
    T_start, T_end,total_neurons_pop1, total_neurons_pop2)

%% shrink the spk index and time by keeping the nonzero ones
spkIndex = spkIndex(spkTime>0)+1;
spkTime = spkTime(spkTime>0);

spkIndex_corr = spkIndex;
spkTime_corr = spkTime;

%%
bin_corr = 1; %ms; changing this is not recommended (need to rescale other numbers throughout)
tEdge = T_start : bin_corr : T_end;
neuronEdge_Pop1 = [index_start_pop1-0.5:1:index_end_pop1+0.5];
neuronEdge_Pop2 = [index_start_pop2-0.5:1:index_end_pop2+0.5];

%% binned spike count for each exc pop
bSpk_Pop1 = histcounts2(spkIndex_corr, spkTime_corr, neuronEdge_Pop1, tEdge);
bSpk_Pop2 = histcounts2(spkIndex_corr, spkTime_corr, neuronEdge_Pop2, tEdge);

%% Find the cross-correlation

tMax = 250;
tMax_bin = tMax/bin_corr;
% subtract the mean for spk train of population 1
sum_spkTrain_Pop1 = sum(bSpk_Pop1,1);
sum_spkTrain_Pop1 = sum_spkTrain_Pop1 - mean(sum_spkTrain_Pop1);

sum_spkTrain_Pop2 = sum(bSpk_Pop2,1);
sum_spkTrain_Pop2 = sum_spkTrain_Pop2 - mean(sum_spkTrain_Pop2);

[total_corr, lags] = xcorr(sum_spkTrain_Pop1, sum_spkTrain_Pop2,tMax_bin); 

%% Normalize the total cross correlation
corrAO_sim = total_corr / ((T_end-T_start)*total_neurons_pop1*total_neurons_pop2);

%% Take the Fourier transform
fmaxHz=100;
T_lim=400;

index_start = find(lags==-100,1);
index_end = find(lags==100,1);

nffsp=floor(2*T_lim*fmaxHz/1000); 
ffsp=(1/(2*T_lim))*[1:nffsp]';
wfsp=ffsp*2*pi;
sfsp=zeros(nffsp,1);
for k=1:nffsp   
   sfsp(k) = trapz(lags(index_start:index_end),...
       exp(-1i*wfsp(k)*lags(index_start:index_end)).*corrAO_sim(index_start:index_end));
end

end

