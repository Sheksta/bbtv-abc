function [neighbours,non_neighbours] = genneighbours(adj)
%% GENNEIGHBOURS
%This function creates a list on neighbours/non-neighbours given an 
%adjacency matrix describing a network.

%Initialising/Pre-allocating outputs


%% Create list of neighbours
[i,j]=find(adj);                                                
neighbours = accumarray(i,j,[size(adj,1), 1],@(x) {sort(x).'});

%% Create list of non-neighbours
adj_NaN=adj;
adj_NaN(logical(eye(size(adj_NaN)))) = 2;
[i,j]=find(adj_NaN==0);                                                
non_neighbours = accumarray(i,j,[size(adj_NaN,1), 1],@(x) {sort(x).'});
end
