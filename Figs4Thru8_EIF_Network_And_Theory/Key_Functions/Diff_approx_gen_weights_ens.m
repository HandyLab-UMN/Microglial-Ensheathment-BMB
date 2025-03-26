%%
% Create the connection matrix used in the diffusion approximation such
% that the out-degrees are fixed
% Note: There are no feedforward connections, since these are approximated
% in the main function
%%
function [wind,wipost,wstr,ensh_realizations,tausyn] = ...
    Diff_approx_gen_weights_ens_v2(Ncells,p0,J,prob_en,s_en,tausyn_default,dt)

%% set up recurrent weight matrix
Ntot = sum(Ncells);
Npop = length(Ncells);

Maxw = round(Ntot*Ntot*0.3); % maximum number of weights in the weight matrix
wind = zeros(Ntot+1,1); % column of w corresponding to the start of the ith neuron's projections                
	
wipost = zeros(Maxw,1);
wstr = zeros(Maxw,1);
ensh_realizations = zeros(Maxw,1);

syncount = 1;
% loop through the populations
for pp = 1:Npop
    
    % loop through each neuron
    for cc=1:Ncells(pp)
        
        starting_Index = sum(Ncells(1:(pp-1)));
        wind(cc + starting_Index) = syncount;
        
        % find which neurons are connected to neuron cc
        for qq = 1:Npop
        
            % probability of a connection
            prob = p0(pp,qq);
            
            % Note: this method fixes the total number of connections to be Ncells(qq)*prob
            iconns = randperm(Ncells(qq),round(Ncells(qq)*prob))+sum(Ncells(1:(qq-1)));

            % Alternative Method: flip weighted coins!            
            % iconns = find(rand(Ncells(qq),1) < prob) + sum(Ncells(1:(qq-1)));
            
            % Find the synapses that are ensheathed
            ensh_realizations(syncount:(syncount+length(iconns)-1)) = ...
                randsample(length(s_en),length(iconns),true,prob_en(pp,:))-1;

            % record the connections
            wipost(syncount:(syncount+length(iconns)-1)) = iconns;
            wstr(syncount:(syncount+length(iconns)-1)) = ...
                J(pp,qq)*(1-s_en(ensh_realizations(syncount:(syncount+length(iconns)-1))+1));
            
            syncount = syncount+length(iconns);
        end
    end
end

% Update the synaptic decay constant
% Synaptic kinetics take the form: 1/tau*exp(-t/tau)*Heaviside(t)
% Indices 1 to length(s_en): Exc constants
% Indices 1+length(s_en) to 2*length(s_en): PV constants
% Indices 1+2*length(s_en) to end: SST constants
tausyn = zeros(3*length(s_en),1);
tausyn_min = tausyn_default*0.4; % i.e., a 60 percent decrease minimum
for ii = 1:length(s_en)
    tausyn(ii) = tausyn_min(1)+(tausyn_default(1)-tausyn_min(1))*(1-s_en(ii));
    tausyn(ii + length(s_en)) = tausyn_min(2)+(tausyn_default(2)-tausyn_min(2))*(1-s_en(ii));
    tausyn(ii + 2*length(s_en)) = tausyn_min(3)+(tausyn_default(3)-tausyn_min(3))*(1-s_en(ii));
end
% repeat these parameters for the surround populations
tausyn = repmat(tausyn,Npop/3,1); 
if max(tausyn<dt)==1
    warning('foo:bar',['Synaptic time constant is smaller than the time step\n' ...
         'Setting those tau_s = dt']);
    tausyn(tausyn<dt)=dt;
end

wind(Ntot+1)= syncount-1;
% reshape these vectors to get rid of unnecessary zeros
wipost = wipost(1:(syncount-1));
wstr = wstr(1:(syncount-1));
ensh_realizations = ensh_realizations(1:(syncount-1));

end