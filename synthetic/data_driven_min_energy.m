function [u_opt,norm_u,err] = data_driven_min_energy(U,Y,sys,T,yf)
% DATA_DRIVEN_MIN_ENRGY - computes data-driven minimum-energy control
%
% Syntax:  [u_opt,norm_u,err] = data_driven_min_energy(U,Y,T,yf)
%
% Inputs:
%    U - input data matrix
%    Y - output data matrix
%    sys - state space system
%    T - control time horizon
%    yf - target state
%
%
% Outputs:
%    u_opt - minimum-energy input sequence
%    norm_u - norm minimum-energy input sequence
%    err_u - error in the final state
%
% Author: Giacomo Baggio
% email: baggio [dot] giacomo [at] gmail [dot] com
% March 2020; Last revision: 01-April-2020

%------------- BEGIN CODE --------------

% number of nodes
n = size(sys.A,1);
% number of inputs
m = size(sys.B,2);
% number of outputs
p = size(sys.C,1);

% compute (output) controllability matrix
C_o = zeros(n,m*T);

C_o(:,1:m) = sys.B;

for k=1:T-1
    
    C_o(:,(m*k+1):(k+1)*m) = sys.A*C_o(:,((k-1)*m+1):m*k);
    
end

C_o = sys.C*C_o;

% data-driven control
K = null(Y);
[W,S,V]=svds(U*K,m*T-p);

% data-driven control
u_opt = U*pinv(Y,1e-8)*yf-U*K*pinv(W*S*V')*U*pinv(Y,1e-8)*yf;

% error on final state
err = norm(yf - C_o*u_opt);
% input norm
norm_u = norm(u_opt);


end

%------------- END OF CODE --------------
