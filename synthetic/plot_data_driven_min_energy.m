% PLOT_DATA_DRIVEN_MIN_ENERGY - Plots error and control energy of minimum-energy 
% model-based and data-driven controls for Erdős–Rényi dynamical networks 
% as a function of the data size
%
% Other m-files required: data_driven_min_energy.m, data_driven_min_energy_approx.m
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

% norm min-energy model-based input
norm_u_opt = zeros(1,length(int_N));
% norm min-energy data-driven input
norm_u_data_driven = zeros(1,length(int_N));
% norm approx min-energy data-driven input
norm_u_data_driven_approx = zeros(1,length(int_N));

% error min-energy model-based input
err_opt = zeros(1,length(int_N));
% error min-energy data-driven input
err_data_driven = zeros(1,length(int_N));
% error approx min-energy data-driven input 
err_data_driven_approx = zeros(1,length(int_N));

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

% compute (output) controllability matrix
C_o = zeros(n,m*T);

C_o(:,1:m) = B;

for k=1:T-1
    
    C_o(:,(m*k+1):(k+1)*m) = sys.A*C_o(:,((k-1)*m+1):m*k);
    
end

C_o = sys.C*C_o;

l=1;

for N = int_N 
    
    err_opt(l) = 0;
    norm_u_opt(l) = 0;
    err_data_driven(l) = 0;
    norm_u_data_driven(l) = 0;
    
    for q=1:Ntrials
       
        % create random data matrices
        U = randn(m*T,N);
        Y = C_o*U;
        
        % data-driven control
        [u_opt_dd,norm_u_dd,err_dd] = data_driven_min_energy(U,Y,sys,T,yf);
        % approx data-driven control
        [u_opt_dd_approx,norm_u_dd_approx,err_dd_approx] = data_driven_min_energy_approx(U,Y,sys,T,yf);
        
        % model-based control
        u_opt = pinv(C_o)*yf;
        
        err_data_driven(l) = err_data_driven(l)+err_dd/Ntrials;
        norm_u_data_driven(l) = norm_u_data_driven(l)+norm_u_dd/Ntrials;
        
        err_data_driven_approx(l) = err_data_driven_approx(l)+err_dd_approx/Ntrials;
        norm_u_data_driven_approx(l) = norm_u_data_driven_approx(l)+norm_u_dd_approx/Ntrials;

        err_opt(l) = err_opt(l)+norm(yf - C_o*u_opt)/Ntrials;
        norm_u_opt(l) = norm_u_opt(l)+norm(u_opt)/Ntrials;
        
    end
    
    l = l+1;
    
end

figure

plot(int_N,log10(err_opt),'-o');
hold on
plot(int_N,log10(err_data_driven),'-*');
hold on
plot(int_N,log10(err_data_driven_approx),'-x');
xlabel('# samples')
ylabel('error in final state')
legend({'optimal','data-driven','data-driven approx'},'Location','best')

figure

plot(int_N,log10(norm_u_opt),'-o');
hold on
plot(int_N,log10(norm_u_data_driven),'-*');
hold on
plot(int_N,log10(norm_u_data_driven_approx),'-x');
xlabel('# samples')
ylabel('control energy')
legend({'optimal','data-driven','data-driven approx'},'Location','best')

%------------- END OF CODE --------------
