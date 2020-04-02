% PLOT_DATA_DRIVEN_VS_SUBSPACE_ID - Compares error of minimum-energy 
% model-based (where the models is estimated via subspace-based methods)
% and data-driven controls for Erdős–Rényi dynamical networks as a 
% function of the network dimension
%
% Other m-files required: data_driven_min_energy.m
%
% Author: Giacomo Baggio
% email: baggio [dot] giacomo [at] gmail [dot] com
% March 2020; Last revision: 01-April-2020

%------------- BEGIN CODE --------------

clc
clear all
close all

int_n = 50:10:200; % network dimension
epsilon = 0.1; % E-R edge density
Ntrials = 50; % number of trials

norm_u_opt = zeros(Ntrials,length(int_n));
norm_u_data_driven = zeros(Ntrials,length(int_n));
norm_u_data_driven_sub = zeros(Ntrials,length(int_n));

err_opt = zeros(Ntrials,length(int_n));
err_data_driven = zeros(Ntrials,length(int_n));
err_data_driven_sub = zeros(Ntrials,length(int_n));

l=1;

for n = int_n
         
    m = n/10; % input dimension
    p = n; % output dimension
    T = 20; % control horizon
    N = m*T+10;  % samples
    yf = randn(p,1); % final state

    for q=1:Ntrials
        
        C_o = zeros(p,m*T);
        A = zeros(n);
        G = graph(A,'omitselfloops');
        [bins,binsize] = conncomp(G);
        
        % ensure that the E-R graph is always connected & controllable
        while min(svd(C_o))<1e-14 && length(binsize)>1
            
            % adjacency matrix E-R graph
            A = randomgraph(n,log(n)/n+epsilon);
            G = graph(A,'omitselfloops');
            [bins,binsize] = conncomp(G);
            A_bin = A;
            A = A/(sqrt(n));
            
            % input matrix
            B = zeros(n,m);
            
            b_temp = randperm(n);
            b = b_temp(1:m);
            
            for qq = 1:length(b)
                
                B(b(qq),qq) = 1;
                
            end
            
            % output matrix
            C = eye(n);
            
            % state space system
            sys = ss(A,B,C,[]);
            
            % compute (output) controllability matrix
            C_o = zeros(n,m*T);
            
            C_o(:,1:m) = B;
            
            for k=1:T-1
                
                C_o(:,(m*k+1):(k+1)*m) = sys.A*C_o(:,((k-1)*m+1):m*k);
                
            end
            
            C_o = sys.C*C_o;
            
        end
        
        % create random data matrices
        U = rand(m*T,N);
        Y = C_o*U;
        
        % data-driven control
        [u_opt_dd,norm_u_dd,err_dd] = data_driven_min_energy(U,Y,sys,T,yf);

        %% model identification via subspace method via controllability matrix
        
        % estimate controllability matrix
        Ctrb_est = Y*pinv(U);
        [Us,Sigma,Vs] = svd(Ctrb_est);
        % select model order
        nr = n;
        
        % estimate C
        C_est = eye(n);
        % estimate B
        B_est = Ctrb_est(:,1:m);
        
        Ctrb_est_p = Ctrb_est(:,m+1:end);
        Ctrb_est_m = Ctrb_est(:,1:(end-m));
        % estimate A
        A_est = Ctrb_est_p*pinv(Ctrb_est_m);
        
        % estimated system
        sys_est = ss(A_est,B_est,C_est,[],-1);

        % compute (output) controllability matrix estimated system
        C_o_est = zeros(n,m*T);
        
        C_o_est(:,1:m) = sys_est.B;
        
        for k=1:T-1
            
            C_o_est(:,(m*k+1):(k+1)*m) = sys_est.A*C_o_est(:,((k-1)*m+1):m*k);
            
        end
        
        C_o_est = sys_est.C*C_o_est;
        
        % model-based control
        u_opt = pinv(C_o_est)*yf;
        
        
        err_data_driven(q,l) = err_dd;
        norm_u_data_driven(q,l) = norm_u_dd;
        
        err_opt(q,l) = norm(yf - C_o*u_opt);
        norm_u_opt(q,l) = norm(u_opt);
        
    end
    
    l = l+1;
     
end

figure

plot(int_n,log10(mean(err_opt,1)),'-o');
hold on
plot(int_n,log10(mean(err_data_driven,1)),'-*');
xlabel('dimension')
ylabel('error in final state')
legend({'optimal','data-driven'},'Location','best')

%------------- END OF CODE --------------
