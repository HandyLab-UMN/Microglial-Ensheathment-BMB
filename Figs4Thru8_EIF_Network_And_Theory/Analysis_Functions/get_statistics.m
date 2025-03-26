%%
% Finds relevant statistics from the spike train s
%%
function [rSim, scc_all,h_all, b_all]=...
    get_statistics(s,T,pInds,Ncells,Tburn)

% Firing rate statistics for all populations
s=s(:,s(1,:)>0);

rSim = zeros(length(Ncells),1);
for ii = 1:length(Ncells)
    rSim(ii)=1000*nnz(s(1,:)>Tburn & s(2,:)>=pInds(ii) & s(2,:)<pInds(ii+1))/(Ncells(ii)*(T-Tburn));
end

%% %%%%%%%%%%%%%%%%%%%
% Estimate spike count correlation distribution
%%%%%%%%%%%%%%%%%%%%%
Ne = Ncells(1);

% counting window size
winsize=250;

% number of cells to sample
nc=4000;

% Edges for histogram
edgest=0:winsize:T;
edgesi=(1:Ne+1)-.01;
edges={edgest,edgesi};

% Find excitatory spikes,
% store into s0, which has the structure
% needed for hist3
Is=find(s(2,:)>0);
s0=s(:,Is)';

% Get 2D histogram of spike indices and times
counts=hist3(s0,'Edges',edges);

% Get rid of edges, 
% the last element in each
% direction is zero
counts=counts(ceil(Tburn/winsize):end-1,1:end-1);

% Find neurons with rates >=2Hz
Igood=find(mean(counts)/winsize>1/1000*0);

% Randomly choose nc cells to use
% their indices are stored in Inds
temp=randperm(numel(Igood),min(nc,numel(Igood)));
Inds=Igood(temp);

% Store their spike counts
counts=counts(:,Inds);

% Compute their pairwise covariances (not correlations; for that use corrcoef
C=cov(counts); 

% Only keep the lower-half of the correlation matrix
% store into SCorree
[II,JJ]=meshgrid(1:size(C,1),1:size(C,1));
SCorree=C(II>JJ);


% Get spike count correlation from all, same pop, and opp pop
scc_all = SCorree;

% Get correlation distributions
[h_all,b_all]=hist(SCorree,100);
h_all=h_all./trapz(b_all,h_all);

end

