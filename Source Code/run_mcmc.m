%%Load data
load('Load data.mat')

%%Covariance Matrix from test runs
cov_rw = [0.0267    0.0002    0.0014   -0.0076   -0.0001   -0.0006;...
    0.0002    0.0758   -0.0293   -0.0007    0.0013   -0.0006;...
    0.0014   -0.0293    0.0160    0.0039   -0.0019    0.0000;...
   -0.0076   -0.0007    0.0039    0.0489   -0.0063    0.0036;...
   -0.0001    0.0013   -0.0019   -0.0063    0.0689   -0.0230;...
   -0.0006   -0.0006    0.0000    0.0036   -0.0230    0.0122];

dist_target = 23; % target tolerance
N = 1e7;  % number of mcmc iterations
theta_init = [0.2283    0.0886    0.0055    0.3294    0.0037    0.0098];  % starting value of chain

%% Run MCMC ABC parameter estimation algorithm
[part_vals, part_s,part_obs,part_sim] = MCMC_bananas_TP(N,boolinfCountQ,border,dist_target,1:37,cov_rw,theta_init);


