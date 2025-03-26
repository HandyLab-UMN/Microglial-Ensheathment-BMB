%%
% Default parameters used in Garcia et al. (2025)
%
% Inputs:
%   surround_diff (scalar): 
%       Difference in orientation between center and surround 
%       Between 0 and 90 degrees
%   contrast (scalar):
%       Intensity of the contrast 
%       Between 0.5 and 1
%   light (cell array):
%       light{1} = 'on' has VIP firing at a low rate
%       light{1} = 'off' has VIP firing at a high rate
%       light = {'range', #} has VIP firing at # (in Hz) 
%   stim_size:
%       Size of stimulus, indicating the number of retinotopic locations
%       'small' (1 center), 'med' (1 center, 1 surround), or 'large' (1
%       center, 2 surrounds)
%
% Written by Gregory Handy, 09/18/2024
%%
function [params] = EIF_params_official_fn(surround_diff,contrast,light,stim_size)

%% Numerical parameters for spiking sims

%Use for accurate approximation of gamma power
% params.T = 1000000; % total time
% params.dt = 0.025;

% Use for faster simulation (to test that system is working)
params.T = 5000; % total time, use this for shorter sims
params.dt = 0.025;
params.NT = round(params.T/params.dt); % total number of time steps
params.recstart = 1000; % time at which to start analysis, to remove transient

%% Population parameters
if strcmp(stim_size,'small') || strcmp(stim_size,'spont')
    params.Npop = 3;
elseif strcmp(stim_size,'med')
    params.Npop = 6;
elseif strcmp(stim_size,'large')
    params.Npop = 9;
    if surround_diff ~=0 
       error('Large stim size only tuned for iso surround condition') 
    end
else
    error('Stim size not recognized')
end
    
% Number of excitatory neurons
params.Ncells = zeros(params.Npop,1); % number of cells in each population: [e, pv, som, vip]
params.Ncells = repmat([4000, 500, 500]',params.Npop/3,1);
params.Ntot = sum(params.Ncells);

params.pinds = cumsum([0;params.Ncells])+1;   % start index of each population (ends with Ntot)
params.rates = zeros(params.Ntot,1);

% Kills the sim if the neuron rate is >50 Hz
% Meaning, every neuron is spiking at >50 Hz
if params.T < 1000000
    params.maxns = round(params.T*params.Ntot*0.1);
else
    % Make sure 20 Hz is a reasonable adjustment!
    params.maxns = round(params.T*params.Ntot*0.025);
end

maxIntSize = 2147483647;
if maxIntSize<params.maxns
   error('The max number of spikes is too large, errors will occurr'); 
end

%% Connectivity parameters
% Center connectivity matrix
p_center = [0.07 0.15 0.10;
            0.05 0.10 0.10;
            0.10 0.00 0.00];  
p_iso = [0.02 0.03 0.08];
p_cross = [0.01 0.005 0.05];
params.p = zeros(params.Npop,params.Npop);
params.p(1:3,1:3) = p_center;
if strcmp(stim_size,'med') || strcmp(stim_size,'large')
    params.p(4:6,4:6) = p_center;
    params.p(4:6,1) = p_iso+(p_cross-p_iso)*surround_diff/90;
    params.p(1:3,4) = params.p(4:6,1);
end
if strcmp(stim_size,'large')
    params.p(7:9,7:9) = p_center;
    params.p(7:9,1) =[0.01 0.015 0.025];
    params.p(7:9,4) =[0.02 0.03 0.05];
    params.p(1:6,7) =[0.01 0.015 0.025 0.02 0.03 0.05];
end
params.p = params.p';

% Connection weights
params.w = 0.48;
params.g = 4;
params.W = params.w*repmat([1, -params.g, -params.g],params.Npop,params.Npop/3);
params.gSyn = params.W';   
params.J = params.gSyn;

%% Neuron-specific parameters

% cell parameters, one value for each population
params.EL = -60*ones(1,params.Npop); %leak voltage, mV
params.vth = 20*ones(1,params.Npop); %threshold voltage, mV
params.vre = -75*ones(1,params.Npop); %reset voltage, mV

params.tau_m = 5.4*ones(1,params.Npop); %Cm./gL, membrane time constant (ms);
params.tau_s = 0.6*ones(1,params.Npop); % synaptic time constant (ms)
params.tau_ref = 1.2*ones(1,params.Npop); %refractory period (ms)
params.DeltaT = 1.4*ones(1,params.Npop); % exponential shape parameter
params.VT = -50*ones(1,params.Npop); % soft threshold
params.vlb = -100; % lower bound for voltage

params.tauDelay = round(repmat([1.8 1.8 1.8],1,params.Npop/3)/params.dt); % synaptic delay; unitless: msec/stepsize
% params.tauDelay = round(repmat([0 0 0],1,params.Npop/3)/params.dt); % synaptic delay; unitless: msec/stepsize
params.tau_d = params.tauDelay*params.dt;
params.maxDelay = max(params.tauDelay);

params.shared_noise_mag = 0.25*ones(1,params.Npop);

%% Mean and standard deviation of feedforward input (background and stimulus)
params.mu_bg = repmat([3 3 7],1,params.Npop/3);
params.sigma_bg = repmat([sqrt(4.5),sqrt(4.5),3],1,params.Npop/3);

params.mu_stim = repmat([3,3,0],1,params.Npop/3)*contrast;
params.sigma_stim = repmat([sqrt(4.5),sqrt(4.5),0],1,params.Npop/3)*sqrt(contrast);

%% VIP parameters

% Number of VIP neurons 
params.N_vip = repmat([0,0,500],1,params.Npop/3);

% total strength of VIP connections
params.J_vip = (-params.g*params.w*2)*params.N_vip;
params.J_sigma_vip = (params.g*params.w*2)^2/(2*params.tau_s(1))*params.N_vip;

% The firing rate of the VIP neuron depending on the light state
if strcmp(light{1},'off')
    r_vip = 4;
elseif strcmp(light{1},'on') % suppresses VIP 
    r_vip = 2;
elseif strcmp(light{1},'range')
    r_vip = light{2};
else
    error('Unclear light state');
end
params.mu_vip(1,:) = params.J_vip.*(r_vip*1e-3);
params.sigma_vip(1,:) = sqrt(params.J_sigma_vip.*(r_vip*1e-3));


%% Corresponding parameters for linear response theory

% Current not used because we only care about w = 0
Tmax = 500;         % Maximum time lag over which to calculate cross-correlations (ms)
dt = 0.5;             % Bin size for which to calculate cross-correlations (ms)
% Generate a vector of frequencies at which to solve for the spectral
% statistics in order to generate cross-correlations with maximum lag Tmax
% and bin size dt.
dw = 1/2/Tmax;
wmax = 1/2/dt;
w = -wmax:dw:(wmax-dw);
ind0 = find(abs(w) < 1e-6);
w(ind0) = 1e-6;

params.ind0 = ind0;
params.omega = w;

params.bins = length(params.omega);

params.dv = 10E-3;         % Membrane potential step to use in solving for the statistics

%%
% in-degrees    
params.I = params.p'.*params.Ncells';
params.J_theory = params.I.*params.W;

end