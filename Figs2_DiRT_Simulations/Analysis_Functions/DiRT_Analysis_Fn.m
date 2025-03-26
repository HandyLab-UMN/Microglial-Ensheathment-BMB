%% 
% Ave the results and fit the curves
%%
function [DiRT_results_ave,alphaFn, JMax,tau_s,rmse,rsquared,...
    pwlMDL_J, pwlMDL_tau] = DiRT_Analysis_Fn(phi_vec,DiRT_results)

% Average across trials
DiRT_results_ave = squeeze(mean(DiRT_results,3));

% We will be fitting these two functions 
alphaFn = @(b,t) b(1)/b(2)^2*t.*exp(-t/b(2));

numPhi = 1:size(phi_vec,2);

JMax = zeros(length(numPhi),1);
tau_s = zeros(length(numPhi),1);
rmse = zeros(length(numPhi),1);
rsquared = zeros(length(numPhi),1);
parfor ii = 1:length(numPhi)

    % Initial guesses [Jmax, tau_s]
    initialGuess = [max(DiRT_results_ave(:,2,ii))*exp(1)*0.1 0.1];

    % Fit the alpha-fn model
    alphaMDL = fitnlm(DiRT_results_ave(:,1,ii),DiRT_results_ave(:,2,ii),...
        alphaFn,initialGuess);

    % Save the values;
    JMax(ii,:) = [alphaMDL.Coefficients.Estimate(1)];
    tau_s(ii,:) = [alphaMDL.Coefficients.Estimate(2)];

    rmse(ii,:) = [alphaMDL.RMSE];
    rsquared(ii,:) = [alphaMDL.Rsquared.Ordinary];
end


sympref('HeavisideAtOrigin', 0);
pwlFn =  @(b,x) b(1)+(b(2).*(x-b(3))).*heaviside(x-b(3));
initialGuess = [JMax(1,1) -JMax(1,1) 0];
pwlMDL_J = fitnlm(phi_vec,JMax(:,1),pwlFn,initialGuess);

notchLoc = pwlMDL_J.Coefficients.Estimate(3);
initialGuess = [tau_s(1,1) -tau_s(1,1)];
pwlFn_v2 =  @(b,x) b(1)+(b(2).*(x-notchLoc)).*heaviside(x-notchLoc);
pwlMDL_tau = fitnlm(phi_vec,tau_s(:,1),pwlFn_v2,initialGuess);


end