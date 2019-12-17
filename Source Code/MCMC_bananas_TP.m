function [part_vals, part_s,part_obs,part_sim] = MCMC_bananas_TP(N,A,adj,dist_target,select_all,cov_rw,theta_init)
%% MCMC_bananas_TP
% Code for performing the MCMC ABC algorithm on the BBTV Forward-Simulation Model. 
% Inputs:
%   
%   N - the number of particles
%   A - The observed data (boolinfQ)
%   adj - adjacency matrix (border)
%   dist_target - the initial discrepancy distance
%   select_all - provisionary input for when only specific months are to be
%   compared. Use select_all = 1:37 for standard simulation.
%   cov_rw - Covariance matrix of the parameters from a previous run.
%   theta_init - Starting value of parameters. Convergence values from
%   prior trial run may be used.

%
% Outputs:
%   part_vals - the particles values
%   part_s - the particles discrepancies
%   part_obs - the observed summary statistics
%   part_sim - the simulated summary statistics
%% 
[neighbours,non_neighbours] = genneighbours(adj);              %generate list of neighbours, and non-neighbours
[~,~,~,~,~,~,part_obs] = smry_new_testing(A,adj,select_all); %generate summary statistic of the observed data

%initialise arrays
part_vals = zeros(N,6);
part_s = zeros(N,1);
part_sim = zeros(N,length(part_obs));

theta_curr = theta_init;
theta_curr_trans = log(-log(theta_curr./(1-theta_curr))); %transpose values to new domain, -inf to inf.
part_sim_curr = zeros(1,length(part_obs));
part_s_curr = 0;

for i = 1:N
    theta_prop_trans = mvnrnd(theta_curr_trans, cov_rw); % generate proposal on transposed domain
    theta_prop = 1./(1 + exp(exp(theta_prop_trans))); % revert proposal to original domain
    
    logprior_curr = sum(-2*log(1 + exp(exp(theta_curr_trans))) + exp(theta_curr_trans) + theta_curr_trans);  %log odds of current distribution
    logprior_prop = sum(-2*log(1 + exp(exp(theta_prop_trans))) + exp(theta_prop_trans) + theta_prop_trans);  %log odds of proposal distribution
    
    if (rand > exp(logprior_prop - logprior_curr )) 
        % early rejection - save stuff and go to next iteration
        part_vals(i,:) = theta_curr;
        part_s(i) = part_s_curr;
        part_sim(i,:) = part_sim_curr;
        continue;
    end
    
    %Simulate proposed distribution
    data_s = simuldata_seasons(A,theta_prop,neighbours,non_neighbours);
    [~,~,~,~,~,~,part_sim_prop] = smry_new_testing(data_s,adj,select_all);
    
    %generate discrepancy measure
    part_summ_prop = (part_obs-part_sim_prop).^2;
    part_s_prop = mean(part_summ_prop,2); 
    
   % accept if proposal within tolerance 
    if (part_s_prop <= dist_target) 
        fprintf('***** ACCEPTED *****\n');
        theta_curr_trans = theta_prop_trans;
        theta_curr = theta_prop;
        part_sim_curr = part_sim_prop;
        part_s_curr = part_s_prop;
    end
    
    part_vals(i,:) = theta_curr;
    part_s(i) = part_s_curr;
    part_sim(i,:) = part_sim_curr;
end

end