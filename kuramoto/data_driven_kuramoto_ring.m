% DATA_DRIVEN_KURAMOTO_RING - compute optimal data-driven control input
% to enforce a prescribed phase-locked pattern in a Kuramoto ring network
%
% Author: Giacomo Baggio
% email: baggio [dot] giacomo [at] gmail [dot] com
% October 2020; Last revision: 16-November-2020

%------------- BEGIN CODE --------------

clear all
close all
clc

%% generate ring network and set Kuramoto parameters

% network size
n = 10;
% control node size
m = 10;
% control nodes
%m_set = [4 5 6];
m_set = [1 2 3 4 5 6 7 8 9 10];
% number of neighbors per node
k = 2;

% define patterns
for i = 1:n

    theta1(i,1) = wrapTo2Pi(0*pi*i/n); % initial phase-locked pattern
    theta2(i,1) = wrapTo2Pi(2*pi*i/n); % final phase-locked pattern    
    
end

% natural frequencies
omega = zeros(1,n)';

% Construct a regular ring lattice: a graph with n nodes, each connected to k neighbors, k/2 on each side
A = zeros(n);
kHalf = k/2;
rows = reshape(repmat([1:n]', 1, k), n*k, 1);
columns = rows+reshape(repmat([[1:kHalf] [n-kHalf:n-1]], n, 1), n*k, 1);
columns = mod(columns-1, n) + 1;
A = sparse(rows, columns, ones(n*k, 1));
A = full(A);

%% generate data

% number of data 
N = 2000;

X0 = [];
X_f = [];
X_bar = [];
U = [];

tspan = [0 0.5]; % control times
h = 0.01; % discretization step
T = tspan(end)/h;

for l=1:N
    
    theta = theta1+0.1*randn(n,1); % initial condition
    
    for t = 1 : tspan(end)/h - 1
        
        u(:,t) = 0.1*randn(m,1);
        
        p=1;
        
        for node = 1 : n
            
            if ismember(node,m_set)
                theta(node,t+1) = theta(node,t) + h*omega(node) + h*u(p,t);
                p=p+1;
            else
                theta(node,t+1) = theta(node,t) + h*omega(node);
            end
            
            for neighbor = 1 : n
                theta(node,t+1) = theta(node,t+1) + h*A(node,neighbor)*sin(theta(neighbor,t) - theta(node,t));
            end
            
        end
    end
    
    X0 = [X0 theta(:,1)];
    U = [U reshape(fliplr(u),[m*(tspan(end)/h-1),1])];
    X_bar = [X_bar reshape(theta(:,2:end-1),[n*(tspan(end)/h-2),1])];
    X_f = [X_f theta(:,end)];
    
end

K_X0 = null(X0,1e-10);
K_U = null(U,1e-10);
xf_c = (theta2-(X_f*K_U*pinv(X0*K_U,1e-10))*theta1);

U = U*K_X0;
X_bar = X_bar*K_X0;
X_f = X_f*K_X0;

% compute data-driven input
K_f = null(X_f,1e-10);
Q = 50*eye(n*(tspan(end)/h-2));
R = 1*eye(m*(tspan(end)/h-1));
L = cholcov(X_bar'*Q*X_bar+U'*R*U);
[W,S,V]=svds(L*K_f,m*(T-1)-n);

u_opt = U*pinv(X_f,1e-10)*xf_c-U*K_f*pinv(W*S*V',1e-10)*L*pinv(X_f,1e-10)*xf_c;

 
%% simulate controlled system
 
u_opt_seq = fliplr(reshape(u_opt,[m,(tspan(end)/h-1)]));

theta = theta1;
tspan = [0 9];

for t = 1 : tspan(end)/h - 1
    N = size(theta,1);
    
    if t<50
        u(:,t) = u_opt_seq(:,t);
    else
        u(:,t) = zeros(m,1);
    end
    
    p=1;
    for node = 1 : n
        
        if ismember(node,m_set) 
            theta(node,t+1) = theta(node,t) + h*omega(node) + h*u(p,t);
            p=p+1;
        else
            theta(node,t+1) = theta(node,t) + h*omega(node);
        end
        
        for neighbor = 1 : n
            theta(node,t+1) = theta(node,t+1) + h*A(node,neighbor)*sin(theta(neighbor,t) - theta(node,t));
        end
        
    end
end


theta_tot = [theta1.*ones(n,100) theta];
figure
plot(1:(tspan(end)/h+100), theta_tot)
xlabel('time');
ylabel('phases');


%------------- END OF CODE --------------
