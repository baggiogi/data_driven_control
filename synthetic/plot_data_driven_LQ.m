% PLOT_DATA_DRIVEN_LQ - Plots error and control energy of optimal model-based 
% and data-driven controls for Erdos-Renyi dynamical networks as a function
% of the data size
%
% Other m-files required: data_driven_LQ.m
%
% Author: Giacomo Baggio
% email: baggio [dot] giacomo [at] gmail [dot] com
% March 2020; Last revision: 01-April-2020

%------------- BEGIN CODE --------------

clc
clear all
close all

n = 100; % network dimension
m = 5; % input dimension
p = 20; % output dimension
T = 10; % control horizon
epsilon = 0.1; % E-R edge density
int_N = 1:100; % samples
Ntrials = 10; % number of trials
yf = randn(p,1); % final state
Q = eye(p*(T-1)); % output cost
R = eye(m*T); % input cost

% norm optimal model-based input
norm_u_opt = zeros(1,length(int_N));
% norm optimal data-driven input
norm_u_data_driven = zeros(1,length(int_N));

% error optimal model-based input
err_opt = zeros(1,length(int_N));
% error optimal data-driveninput
err_data_driven = zeros(1,length(int_N));

% generate E-R nework
A = randomgraph(n,log(n)/n+epsilon);
A_bin = A;
A = A/sqrt(n);

% generate input matrix
B = zeros(n,m);

b_temp = randperm(n);
b = b_temp(1:m);

for q = 1:length(b)
    
    B(b(q),q) = 1;
    
end

% generate output matrix
C = zeros(p,n);

c_temp = randperm(n);
c = c_temp(1:p);

for q = 1:length(c)
    
    C(q,c(q)) = 1;
    
end

% state space system
sys = ss(A,B,C,[]);

% compute T-steps system Hankel matrix
H = zeros(p*T,m*T); 

for r=1:T
    for k=1:T
        if k > T-r
            H((r-1)*p+1:r*p,(k-1)*m+1:k*m) = sys.C*sys.A^(r-T+k-1)*sys.B;
        end
    end
end

H_bar = H(1:p*(T-1),:);

% last block row Hankel matrix
H_f = H(p*(T-1)+1:end,:);

l=1;

for N = int_N 
    
    err_opt(l) = 0;
    norm_u_opt(l) = 0;
    err_data_driven(l) = 0;
    norm_u_data_driven(l) = 0;
    
    for q=1:Ntrials
       
        % create random data matrices
        U = randn(m*T,N);
        Y = H*U;
        Y_bar = Y(1:(T-1)*p,:);
        Y_f = Y((T-1)*p+1:end,:);
        
        % data-driven control
        [u_opt_dd,norm_u_dd,err_dd] = data_driven_LQ(U,Y,sys,T,yf,Q,R);
        
        % model-based control
        K_Hf = null(H_f);
        L_m = cholcov(H_bar'*Q*H_bar+R);
        u_opt = pinv(H_f,1e-8)*yf-K_Hf*pinv(L_m*K_Hf,1e-8)*L_m*pinv(H_f,1e-8)*yf;
        
        err_data_driven(l) = err_data_driven(l)+err_dd/Ntrials;
        norm_u_data_driven(l) = norm_u_data_driven(l)+norm_u_dd/Ntrials;

        err_opt(l) = err_opt(l)+norm(yf - H_f*u_opt)/Ntrials;
        norm_u_opt(l) = norm_u_opt(l)+norm(u_opt)/Ntrials;
        
    end
    
    l = l+1;
    
end

figure

plot(int_N,log10(err_opt),'-o');
hold on
plot(int_N,log10(err_data_driven),'-*');
xlabel('# samples')
ylabel('error in final state')
legend({'optimal','data-driven'},'Location','best')

figure

plot(int_N,log10(norm_u_opt),'-o');
hold on
plot(int_N,log10(norm_u_data_driven),'-*');
xlabel('# samples')
ylabel('control energy')
legend({'optimal','data-driven'},'Location','best')

%------------- END OF CODE --------------
