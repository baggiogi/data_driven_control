function final = smallworld(n, k, beta)
% SMALLWORLD - Generate a small world graph using the "Watts and Strogatz model" as
% described in Watts, D.J.; Strogatz, S.H.: "Collective dynamics of
% 'small-world' networks."
% A graph with n*k/2 edges is constructed, i.e. the nodal degree is n*k for
% every node.
%
% Inputs:
% n: number of nodes of the graph to be generated
% k: mean degree (assumed to be an even integer)
% beta: rewiring probability
%
% Outputs:
% final: sparse symmetric adjacency matrix representing the generated graph

%------------- BEGIN CODE --------------

% Construct a regular lattice: a graph with n nodes, each connected to k
% neighbors, k/2 on each side.
kHalf = k/2;
rows = reshape(repmat([1:n]', 1, k), n*k, 1);
columns = rows+reshape(repmat([[1:kHalf] [n-kHalf:n-1]], n, 1), n*k, 1);
columns = mod(columns-1, n) + 1;
B = sparse(rows, columns, ones(n*k, 1));
A = sparse([], [], [], n, n);

% With probability beta rewire an edge avoiding loops and link duplication.
% Until step i, only the columns 1:i are generated making implicit use of A's
% symmetry.
for i = [1:n]
    % The i-th column is stored full for fast access inside the following loop.
    col= [full(A(i, 1:i-1))'; full(B(i:end, i))];
    for j = i+find(col(i+1:end))'
        if (rand()<beta)
            col(j)=0;
            k = randi(n);
            while k==i || col(k)==1
                k = randi(n);
            end
            col(k) = 1;
        end
    end
    A(:,i) = col;
end

% A is not yet symmetric: to speed things up, an edge connecting i and j, i < j
% implies A(i,j)==1, A(i,j) might be zero.
T = triu(A);
final = T+T';

end 

%------------- END OF CODE --------------

