function [iterations_array] = simuldata_seasons(A,theta_define,neighbours,non_neighbours)
%% SIMULDATA
%This function simulates an infection spreading in a farm, given
%parameters.
%% Inputs:
%A = Observed data [38 * 54 Array]
%adj = Adjacency Matrix [ 54 * 54 Array]
%thetas = Particle Values
%iterate_simul  = Number of simulations
% Outputs:
%iterations_array = simulated data

%% Initialising/Pre-allocating outputs

%Reshape theta_define to be individual parameters along rows, and columns
theta_define = reshape(theta_define, [3,2]);

n = size(A,1);           %    # of time steps
side = size(A,2);       %     # of subsections

iterations_array = zeros(n,side);
iterations_array(1,:,:) = A(1,:);          %All simulations will start with initial config. from Obs data - first time step.

%% Simulation   
    for t=1:n-1
        season = mod(t,12);
        season(season==0)=12;
        if (season >=1 && season <=3 ) || (season >=10 && season <=12 ) 
            season_opt = 1;
        else 
            season_opt = 2;
        end
        
        
        %infected cells at time t
        idx = find(iterations_array(t,:));
        
        %create array for cells with an infected neighbour
        idx_neighbours_inf = [neighbours{idx}];
        idx_neighbours = idx_neighbours_inf;
        idx_neighbours(ismember(idx_neighbours_inf,idx)) = [];
        
        
        %create array for uninfected cells that are not neighbouring at least one
        %infected cell
        idx_non_neighbours_inf = [non_neighbours{idx}];
        idx_non_neighbours = idx_non_neighbours_inf;
        idx_non_neighbours(ismember(idx_non_neighbours_inf,idx)) = [];
        
        for i = 1:length(idx)
            %recovery step - theta 0
            iterations_array(t+1,idx(i)) = rand < 1-theta_define(1,season_opt);
        end
        
        for i = 1:length(idx_neighbours)
            %local infection - theta 1   
            iterations_array(t+1,idx_neighbours(i)) = max(iterations_array(t+1,idx_neighbours(i)) , rand < theta_define(2,season_opt));
        end
        
        for i = 1:length(idx_non_neighbours)
            %ranged infection - theta 2
            iterations_array(t+1,idx_non_neighbours(i)) = max(iterations_array(t+1,idx_non_neighbours(i)) , rand < theta_define(3,season_opt));
        end
        
        
    end
    
end
