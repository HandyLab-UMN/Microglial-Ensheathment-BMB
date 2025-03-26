%%
% Default parameters values (Table 1 in Handy and Borisyuk, 2023)
% All values are saved in the structure DiRT_params to be used elsewhere
%   \phi: fraction of the cleft blocked by protruding astrocyte
%%
function DiRT_params = params_DiRT_2Dcleft(phi,N_rec,tau_r)

%% Key simulation parameters that can easily be adjusted 
% Particle details
DiRT_params.N_NT = 1000;  % number of neurotransmitters released
DiRT_params.D = 1;        % diffusion coefficient

%% Locations of the boundaries

% Adjust for phi
cleft_x_total_default = 1.0;
if phi<-2 || phi >= 1
   error('Parameter phi needs to be in the domain [-2 1). Current value: %0.3f', phi)  
end
% When astro protrusion into the cleft occurs (i.e., phi >0), both the 
% cleft and domain boundaries have to be adjusted due to logic in the mex code
DiRT_params.cleft_x_total = cleft_x_total_default*(1-phi*(phi>0));
DiRT_params.domain_x_total = cleft_x_total_default*(1-phi);

% Dimensions of the cleft
DiRT_params.cleft_x_min = -DiRT_params.cleft_x_total/2;
DiRT_params.cleft_x_max = DiRT_params.cleft_x_total/2;

DiRT_params.cleft_y_total = 0.1;
DiRT_params.cleft_y_min = -DiRT_params.cleft_y_total/2;
DiRT_params.cleft_y_max = DiRT_params.cleft_y_total/2;

% Total dimensions of the domain
DiRT_params.domain_x_min = -DiRT_params.domain_x_total/2;
DiRT_params.domain_x_max = DiRT_params.domain_x_total/2;

DiRT_params.domain_y_total = 10;
DiRT_params.domain_y_min = -DiRT_params.domain_y_total/2;
DiRT_params.domain_y_max = DiRT_params.domain_y_total/2;

%% Capture region locations

DiRT_params.N_rec = N_rec; % number of postsynaptic receptors 

% Total area spanned by the receptors (postsynaptic density)
DiRT_params.capture_region_length = 0.5;
DiRT_params.capture_region_min = -DiRT_params.capture_region_length/2;
DiRT_params.capture_region_max = DiRT_params.capture_region_length/2;

% length of individual receptors
DiRT_params.rec_individual=(DiRT_params.capture_region_max-DiRT_params.capture_region_min)/DiRT_params.N_rec;
for j = 1:DiRT_params.N_rec
    DiRT_params.receptors(j,1) = DiRT_params.capture_region_min+DiRT_params.rec_individual*(j-1);
    DiRT_params.receptors(j,2) = DiRT_params.capture_region_min+DiRT_params.rec_individual*(j);
end


%%

% Receptor details
DiRT_params.tau_r = tau_r;  % recharge time (i.e., recharge rate is 1/tau_r)

% Numerical details
DiRT_params.num_trials = 10;
DiRT_params.max_time = 4; % max time to run the simulation,

% Define where the particles should start
% Gaussian with center at (x,y,z) and standard deviation 0.01
DiRT_params.x_start_loc = 0;
DiRT_params.y_start_loc = DiRT_params.cleft_y_max;


%% No need to adjust the parameters below this point
% Note: simulation will automatically end once all particles have been
% removed and receptors have returned to the open state
DiRT_params.dt = 0.00001; % time step for diffusion update

DiRT_params.dt_saved = 0.0001;
DiRT_params.print_num = floor(DiRT_params.dt_saved/DiRT_params.dt); % when to save the state of the system
% DiRT_params.max_time_points = DiRT_params.print_num*DiRT_params.max_time;
DiRT_params.max_time_points = round(DiRT_params.max_time/DiRT_params.dt_saved);

DiRT_params.t = DiRT_params.dt_saved:DiRT_params.dt_saved:DiRT_params.max_time;

% affinity for the receptor;
DiRT_params.affinity = sqrt(pi)*sqrt(DiRT_params.dt);

end