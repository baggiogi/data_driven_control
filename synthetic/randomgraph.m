function  final= randomgraph(n,p)
% RANDOMGRAPH - generates an undirected random graph (without self-loops) of size n (as
% described in the Erdoes-Renyi model)
%
% Inputs:
%    n - number of nodes
%    p - probability that node i and node j, i != j, are connected by an edge
%
% Outputs:
%    final - nxn full symmetric adjacency matrix representing the generated graph
%

%------------- BEGIN CODE --------------

G = rand(n,n) < p;
G = triu(G,1);
G = G + G';
final = G;

end 

%------------- END OF CODE --------------
