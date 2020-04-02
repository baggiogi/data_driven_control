% PLOT_DATA_COMPUTE_TIME_DATA_DRIVEN - Plots the computational times and errors
% when computing minimum-energy controls using the data-driven strategies and 
% model-based ones for Erdős–Rényi dynamical networks as a function of 
% the network dimension
%
% Author: Giacomo Baggio
% email: baggio [dot] giacomo [at] gmail [dot] com
% March 2020; Last revision: 01-April-2020

%------------- BEGIN CODE --------------

clc
clear all
close all

int_n = 1000:1000:10000; % network dimension
epsilon = 0.05; % E-R edge density

% computational times
elapsedTime_mb = zeros(1,length(int_n));
elapsedTime_dd = zeros(1,length(int_n));
elapsedTime_dd_approx = zeros(1,length(int_n));

% errors
err_mb = zeros(1,length(int_n));
err_dd = zeros(1,length(int_n));
err_dd_approx = zeros(1,length(int_n));

l = 1;

for n = int_n
   
    n
    m = 0.01*n;  % input dimension
    p = 0.01*n; % output dimension
    T = 0.01*n + 10; % control horizon
    N = m*T; % samples
    yf = randn(p,1); % target state
    
    % adjacency matrix E-R graph
    A = randomgraph(n,log(n)/n+epsilon);
    A = A/(norm(A)+0.01);
    
    % input matrix
    B = zeros(n,m);
    
    b_temp = randperm(n);
    b = b_temp(1:m);
    
    for q = 1:length(b)
        
        B(b(q),q) = 1;
        
    end
    
    % output matrix
    C = zeros(p,n);
    
    c_temp = randperm(n);
    c = c_temp(1:p);
    
    for q = 1:length(c)
        
        C(q,c(q)) = 1;
        
    end 
        
    %% model-based computation
    
    % state space system
    sys = ss(A,B,C,[]);

    tic 
    % compute (output) controllability matrix
    C_o = zeros(n,m*T);
    
    C_o(:,1:m) = B;
    
    for k=1:T-1
        
        C_o(:,(m*k+1):(k+1)*m) = sys.A*C_o(:,((k-1)*m+1):m*k);
        
    end
    
    C_o = sys.C*C_o;
    
    % minimum-energy model-based control
    u = pinv(C_o)*yf;
    
    elapsedTime_mb(l) = toc;
    
    err_mb(l) = norm(yf - C_o*u);
    
    %% data-driven computation
    
    U = randn(m*T,p);
    Y = C_o*U;

    tic
    
    % minimum-energy data-driven control
    u_dd = pinv(Y*pinv(U))*yf;
    
    elapsedTime_dd(l) = toc;
    
    err_dd(l) = norm(yf - C_o*u_dd);
    
    tic
    
    % approximate minimum-energy data-driven control
    u_dd_approx = U*pinv(Y)*yf;
    
    elapsedTime_dd_approx(l) = toc;
    
    err_dd_approx(l) = norm(yf - C_o*u_dd_approx);
    
    l = l+1;
    
end

figure

plot(int_n,log10(elapsedTime_mb));
hold on
plot(int_n,log10(elapsedTime_dd));
hold on
plot(int_n,log10(elapsedTime_dd_approx));
xlabel('dimension') 
ylabel('computational time (sec)') 
legend({'optimal','data-driven','data-driven approx'},'Location','best')

figure

plot(int_n,log10(err_mb));
hold on
plot(int_n,log10(err_dd));
hold on
plot(int_n,log10(err_dd_approx));
xlabel('dimension')
ylabel('error in final state')
legend({'optimal','data-driven','data-driven approx'},'Location','best')
