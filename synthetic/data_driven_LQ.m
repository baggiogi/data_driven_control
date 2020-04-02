function [u_opt,norm_u,err] = data_driven_LQ(U,Y,sys,T,yf,Q,R)
% DATA_DRIVEN_LQ - computes data-driven optimal control with linear
% quadratic cost function
%
% Syntax:  [u_opt,norm_u,err] = data_driven_LQ(U,Y,sys,T,yf,Q,R)
%
% Inputs:
%    U - input data matrix
%    Y - (complete) output data matrix
%    sys - state space system
%    T - control time horizon
%    yf - target state
%    Q - output cost
%    R - input cost
%
% Outputs:
%    u_opt - optimal input sequence
%    norm_u - norm optimal input sequence
%    err_u - error in the final state
%
% Author: Giacomo Baggio
% email: baggio [dot] giacomo [at] gmail [dot] com
% March 2020; Last revision: 01-April-2020

%------------- BEGIN CODE --------------

% number of outputs
p = size(sys.C,1);
% number of inputs
m = size(sys.B,2);

% compute T-steps system Hankel matrix
H = zeros(p*T,m*T); 

for r=1:T
    for k=1:T
        if k > T-r
            H((r-1)*p+1:r*p,(k-1)*m+1:k*m) = sys.C*sys.A^(r-T+k-1)*sys.B;
        end
    end
end

% last block row Hankel matrix
H_f = H(p*(T-1)+1:end,:);

% output data in the interval [1,T-1]
Y_bar = Y(1:(T-1)*p,:);
% output data at time T
Y_f = Y((T-1)*p+1:end,:);

K_f = null(Y_f);
L = cholcov(Y_bar'*Q*Y_bar+U'*R*U);
[W,S,V]=svds(L*K_f,m*T-p);

% data-driven control
u_opt = U*pinv(Y_f,1e-8)*yf-U*K_f*pinv(W*S*V')*L*pinv(Y_f,1e-8)*yf;

% error on final state
err = norm(yf - H_f*u_opt);
% input norm
norm_u = norm(u_opt);

end

%------------- END OF CODE --------------
