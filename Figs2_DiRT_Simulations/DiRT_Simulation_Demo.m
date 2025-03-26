%%
% This code runs DiRT simulations
%   This code interfaces with mex for speed
%   Also makes use of parfor across trials (see inner loop)
%   Runs over specified values of phi (i.e., fraction of the cleft blocked 
%       by protruding astrocyte)
%%
clear; close all; clc;

% Add that folder plus all subfolders to the path.
restoredefaultpath;
folder = fileparts(which('DiRT_Simulation_Demo.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Compile the mex code (only needs to happen once)
fprintf('Compiling mex code \n');
mex ./Mex_Code/DiRT_NeuroAstro_mex.c

%% Specified the values of phi for the simulations
% phi_vec: fraction of the cleft blocked by protruding astrocyte
phi_vec = -1:0.1:0.9; % full vector for phi

N_rec = 50; % Number of avaliable receptors
tau_r = 0.1; % receptor recharge time

filename=sprintf('DiRT_Results_N%d_tau%.2f_v2.mat',N_rec,tau_r);

%% Preallocate space
paramsDefault = params_DiRT_2Dcleft(0,N_rec,tau_r);
DiRT_results = zeros(length(paramsDefault.t)+1,2,paramsDefault.num_trials,length(phi_vec));

%% Loop over all phi values
for ii = 1:length(phi_vec)
    
    % Load the parameters for this trial
    params = params_DiRT_2Dcleft(phi_vec(ii),N_rec,tau_r);

    % Seed for the random number generator in the mex code
    randSeed = 833141*ii;
    
    %% Loop over trials 
    % If unable to run in parallel, simply change this to a for-loop
    str1 = sprintf(' %c = %.3f', 981,phi_vec(ii));
    fprintf(strcat('Running the DiRT trials for',str1));
    fprintf(' (%d out of %d values)\n',ii, length(phi_vec));
    output = zeros(params.max_time_points,4,params.num_trials);
    %%
    tic;
    parfor hh = 1:params.num_trials
        % Run the simulation
        [output(:,:,hh)] = DiRT_NeuroAstro_mex(params.max_time,params.dt, params.D, ...
            params.domain_x_min,params.domain_x_max,params.cleft_x_min,params.cleft_x_max,...
            params.domain_y_min,params.domain_y_max,params.cleft_y_min,params.cleft_y_max,...
            params.capture_region_min,params.capture_region_max,...
            params.x_start_loc,params.y_start_loc, params.N_NT, params.tau_r, ...
            params.N_rec, params.receptors, params.print_num, params.max_time_points,...
            params.affinity,randSeed+hh);
    end
    toc;
    
    % Clean up the output a little (can't be in the parfor loop)
    % Some trials might end before max_time is reached (i.e. all particles have 
    % escaped and all receptors recharged). This repeats the end state for the
    % remaining time bins
    for hh = 1:params.num_trials
        temp_index = find(output(:,1,hh)==0,1);
        if ~isempty(temp_index)
            remaining_slots = length(params.t)-temp_index+1;
            output(temp_index:end,2:end,hh) = ...
                repmat(output(temp_index-1,2:end,hh),remaining_slots,1);
        end
        
        % fills in the time column (completely)
        output(:,1,hh) = params.t;
    end
    
    % Save the core results (time and num. of act. receptors)
    % Also concatenate the time = 0 time point
    DiRT_results(:,1,:,ii) = [zeros(1, 1,params.num_trials); output(:,1,:)];
    DiRT_results(:,2,:,ii) = [zeros(1, 1,params.num_trials); params.N_rec - output(:,4,:)];
    
end

%% Save the results
save(strcat('./DiRT_Sim_Data/',filename),'DiRT_results','params','phi_vec','-v7.3')

%% Analyze the output

analyzeDirt_alphaFn
