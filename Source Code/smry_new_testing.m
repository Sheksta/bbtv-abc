function [S, S10,S01_NOT, S01_REDUNDANT,INF,S01,TOTAL] = smry_new_testing(A,adj,select_all)
%Creates summary statistics from provided data,
%given adjacency matrix,data,select_all
for t = 1:37
    season = mod(t,12);
    season(season==0)=12;
    if (season >=1 && season <=3 ) || (season >=10 && season <=12 )
        season_opt(t) = [1];
    else
        season_opt(t) = [0];
    end
end


%Proportion of infected nodes <-- Summary Statistic
    results = sum(A,2);
    INF = (sum(results));
    S = (results([select_all,select_all(end)+1]).'/INF).*100;

%number of infections;
%Vector of length(t)-1 that lists the proportion of cells that were
%infected at time t, and then recovered
    results = sum(diff(A,1,1)==-1,2);
    S10 = results(select_all).';
    S10_sum = sum(S10(season_opt==1));
    S10_wint = sum(S10(season_opt==0));
%     /total_prop;
%Vector of length(t)-1 that lists the proportion of cells that were
%not infected at time t, but then got infected
    results = sum(diff(A,1,1)==1,2);    
    S01 = results(select_all).';
    S01_sum = sum(S01(season_opt==1));
    S01_wint = sum(S01(season_opt==0));
%     S01 = sum(results(select_all),1);
%     /total_prop;

%S01_NOT_obs: Vector of length(t)-1 that lists the proportion of ISOLATED cells that 
%were not infected at time t, but then got infected

    % functions that represent the conditions
    cond1 = @(t,i) A(t,i)==0 ... 
        && A(t+1,i) == 1;     % true if node i is 0 at t, and 1 at t+1
    cond2 = @(t,i) all(adj(i, A(t,:)==1)==0) ; % true if all nodes 1 are not neighbours of node i ; % true if all nodes 1 are not neighbours of node i
    Ni = size(adj,1) ;     % number of nodes
    Nt = size(A,1) - 1 ;   % number of times to check (-1!)
    % engine
    [tt, ii] = ndgrid(1:Nt, 1:Ni) ; 
    R = arrayfun(@(t,i) cond1(t,i) && cond2(t,i) , tt, ii);
    results = sum(R,2);
    S01_NOT = results(select_all).';
    S01_NOT_sum = sum(S01_NOT(season_opt==1));
    S01_NOT_wint = sum(S01_NOT(season_opt==0));
%     S01_NOT = sum(results(select_all),1).';
%     /total_prop.';

%Vector of length(t)-1 that lists the proportion of NON-ISOLATED cells that
%were not infected at time t, but then got infected
   S01_REDUNDANT=S01-S01_NOT;
   S01_RDNT_sum = sum(S01_REDUNDANT(season_opt==1));
   S01_RDNT_wint = sum(S01_REDUNDANT(season_opt==0));    
   TOTAL = [S,S10_sum,S10_wint,S01_NOT_sum,S01_NOT_wint,S01_RDNT_sum,S01_RDNT_wint,S01_sum,S01_wint,INF];
end