% DATA_DRIVEN_FAULT_RECOVERY - compute optimal data-driven control input
% for fault recovery in the New England 10-unit 39-bus Test case
%
% Other m-files required: NE_test_parameters.m
%
% Author: Giacomo Baggio
% email: baggio [dot] giacomo [at] gmail [dot] com
% July 2020; Last revision: 26-July-2020

%------------- BEGIN CODE --------------

clc 
clear all
close all

global H P D E Y3 Y32 T_fault_in T_fault_end

% load grid parameters
NE_test_parameters

% initial stable state
y0 = [0.1564 0.1806 0.1631 0.3135 0.1823 0.1849 0.1652 0.2953 -0.06165 0 0 0 0 0 0 0 0 0]';
y0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

% sampling time
delta_t = 2.5e-4;
% final simulation time
T_sim = 15;
% vector of simulation times
tspan = 0:delta_t:T_sim;
% fault initial time
T_fault_in = 2;
% fault final time
T_fault_end = 2.525;
% number of inputs
m = 9;
% number of states
n = 18;
% control horizon
T = 0.1;
% final stable state
xf = y0;
% control starting time
t_c = 3;

k = 1;

%% simulate fault (discretized swing equations)

for t = tspan
    yy(:,k) = swing(t,delta_t,y0,zeros(1,m));
    y0 = yy(:,k);
    k = k+1;
end

% initial state control
x0 = yy(:,round(t_c/delta_t));

% plot results
figure(1)
subplot(2,2,1)
plot(tspan,yy(1:9,:))
hold on
fill([T_fault_in T_fault_in T_fault_end T_fault_end], [-20 200 200 -20], 'r','FaceAlpha',0.25,'LineStyle','none');
ylabel('phase');
xlabel('time t');
legend('gen 10','gen 2','gen 3','gen 4','gen 5','gen 6','gen 7','gen 8','gen 9')
subplot(2,2,3)
plot(tspan,yy(10:18,:))
hold on
fill([T_fault_in T_fault_in T_fault_end T_fault_end], [-20 60 60 -20], 'r','FaceAlpha',0.25,'LineStyle','none');
ylabel('frequency');
xlabel('time t');
legend('gen 10','gen 2','gen 3','gen 4','gen 5','gen 6','gen 7','gen 8','gen 9')


%% generate data for control

% number of samples
N = 5000;
% variance of inital states
sigma = 0.01;
% data matrices
U = [];
X0 = [];
X1T = [];
XT = [];

for i=1:N
    
    k = 1;
    
    % random input
    u = 1e-3*randn(m,round(T/delta_t));
    % random initial state
    y0 = [0.1564+sigma*randn 0.1806+sigma*randn 0.1631+sigma*randn 0.3135+sigma*randn 0.1823+sigma*randn 0.1849+sigma*randn 0.1652+sigma*randn 0.2953+sigma*randn -0.06165+sigma*randn 0+sigma*randn 0+sigma*randn 0+sigma*randn 0+sigma*randn 0+sigma*randn 0+sigma*randn 0+sigma*randn 0+sigma*randn 0+sigma*randn]';
    
    y0 = [sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn sigma*randn]';
    y0_data = y0;
    
    for t = 0:delta_t:T-delta_t
        y_data(:,k) = swing(t,delta_t,y0_data,u(:,k));
        y0_data = y_data(:,k);
        k = k+1;
    end
    
    U = [U reshape(fliplr(u),[m*round(T/delta_t),1])];
    X0 = [X0 y0];
    X1T = [X1T reshape(fliplr(y_data(:,1:round(T/delta_t-1))),[n*round(T/delta_t-1),1])];
    XT = [XT y_data(:,end)];
end


%% compute data-driven control

U = U(:,1:3750);
XT = XT(:,1:3750);
X0 = X0(:,1:3750);
X1T = X1T(:,1:3750);

K_X0 = null(X0,1e-10);
K_U = null(U,1e-10);
xf_c = (xf-(XT*K_U*pinv(X0*K_U,1e-10))*x0);

U = U*K_X0;
X1T = X1T*K_X0;
XT = XT*K_X0;

K_U = null(U,1e-10);
K_XT = null(XT,1e-10);
Q = 5e-3;
R = 1;

L2 = Q*((X1T)'*(X1T))+R*(U'*U);
L = cholcov(L2);
[W,S,V]=svds(L*K_XT,m*round(T/delta_t)-n);

u_opt = U*pinv(XT,1e-10)*xf_c-U*K_XT*pinv(W*S*V',1e-10)*L*pinv(XT,1e-10)*xf_c;
u_opt_seq = fliplr(reshape(u_opt,[m,round(T/delta_t)]));

%% apply control for fault recovery

u_opt_seq = [zeros(9,round(t_c/delta_t+1)) u_opt_seq zeros(9,round(500/delta_t))];
y0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

k = 1;
tspan = 0:delta_t:500;

for t = tspan
    yy(:,k) = swing(t,delta_t,y0,u_opt_seq(:,k));
    y0 = yy(:,k);
    k = k+1;
end

% plot results
figure(1)
subplot(2,2,2)
plot(0:delta_t:T_sim,yy(1:9,1:round(T_sim/delta_t)+1))
hold on
fill([T_fault_in T_fault_in T_fault_end T_fault_end], [-10 15 15 -10], 'r','FaceAlpha',0.25,'LineStyle','none');
hold on
fill([t_c t_c t_c+T t_c+T], [-10 15 15 -10], 'g','FaceAlpha',0.25,'LineStyle','none');
ylabel('phase');
xlabel('time t');
legend('gen 10','gen 2','gen 3','gen 4','gen 5','gen 6','gen 7','gen 8','gen 9')
subplot(2,2,4)
plot(0:delta_t:T_sim,yy(10:18,1:round(T_sim/delta_t)+1))
hold on
fill([T_fault_in T_fault_in T_fault_end T_fault_end], [-10 15 15 -10], 'r','FaceAlpha',0.25,'LineStyle','none');
hold on
fill([t_c t_c t_c+T t_c+T], [-10 15 15 -10], 'g','FaceAlpha',0.25,'LineStyle','none');
ylabel('frequency');
xlabel('time t');
legend('gen 10','gen 2','gen 3','gen 4','gen 5','gen 6','gen 7','gen 8','gen 9')

%asymptotic behavior
figure(2)
subplot(2,1,1)
plot(tspan,yy(1:9,:))
hold on
fill([2 2 2.4 2.4], [-10 15 15 -10], 'r','FaceAlpha',0.25,'LineStyle','none');
hold on
fill([t_c t_c t_c+T t_c+T], [-10 15 15 -10], 'g','FaceAlpha',0.25,'LineStyle','none');
ylabel('phase');
xlabel('time t');
legend('gen 10','gen 2','gen 3','gen 4','gen 5','gen 6','gen 7','gen 8','gen 9')
subplot(2,1,2)
plot(tspan,yy(10:18,:))
hold on
fill([2 2 2.4 2.4], [-10 15 15 -10], 'r','FaceAlpha',0.25,'LineStyle','none');
hold on
fill([t_c t_c t_c+T t_c+T], [-10 15 15 -10], 'g','FaceAlpha',0.25,'LineStyle','none');
ylabel('frequency');
xlabel('time t');
legend('gen 10','gen 2','gen 3','gen 4','gen 5','gen 6','gen 7','gen 8','gen 9')


function yp= swing(t,delta_t,y,u)
% simulate power generators swing dynamics
    
    global H P D E Y3 Y32 T_fault_in T_fault_end
    
    yp = zeros(18,1);
    
     Y = Y3;
     k = [delta_t*(1/H(2,2))*(P(2)-real(Y(2,2))*E(2)^2-E(2)*E(1)*real(Y(2,1))-E(2)*E(3)*real(Y(2,3))- E(2)*E(4)*(real(Y(2,4)))-E(2)*E(5)*(real(Y(2,5)))-E(2)*E(6)*(real(Y(2,6)))-E(2)*E(7)*(real(Y(2,7)))-E(2)*E(8)*(real(Y(2,8)))-E(2)*E(9)*(real(Y(2,9)))-E(2)*E(10)*(real(Y(2,10))));
        delta_t*(1/H(3,3))*(P(3)-real(Y(3,3))*E(3)^2-E(3)*E(1)*(real(Y(3,1)))-E(3)*E(2)*(real(Y(3,2)))-E(3)*E(4)*(real(Y(3,4)))-E(3)*E(5)*(real(Y(3,5)))-E(3)*E(6)*(real(Y(3,6)))-E(3)*E(7)*(real(Y(3,7)))-E(3)*E(8)*(real(Y(3,8)))-E(3)*E(9)*(real(Y(3,9)))-E(3)*E(10)*(real(Y(3,10))));
        delta_t*(1/H(4,4))*(P(4)-real(Y(4,4))*E(4)^2-E(4)*E(1)*(real(Y(4,1)))-E(4)*E(2)*(real(Y(4,2)))-E(4)*E(3)*(real(Y(4,3)))-E(4)*E(5)*(real(Y(4,5)))-E(4)*E(6)*(real(Y(4,6)))-E(4)*E(7)*(real(Y(4,7)))-E(4)*E(8)*(real(Y(4,8)))-E(4)*E(9)*(real(Y(4,9)))-E(4)*E(10)*(real(Y(4,10))));
        delta_t*(1/H(5,5))*(P(5)-real(Y(5,5))*E(5)^2-E(5)*E(1)*(real(Y(5,1)))-E(5)*E(2)*(real(Y(5,2)))-E(5)*E(3)*(real(Y(5,3)))-E(5)*E(4)*(real(Y(5,4)))-E(5)*E(6)*(real(Y(5,6)))-E(5)*E(7)*(real(Y(5,7)))-E(5)*E(8)*(real(Y(5,8)))-E(5)*E(9)*(real(Y(5,9)))-E(5)*E(10)*(real(Y(5,10))));
        delta_t*(1/H(6,6))*(P(6)-real(Y(6,6))*E(6)^2-E(6)*E(1)*(real(Y(6,1)))-E(6)*E(2)*(real(Y(6,2)))-E(6)*E(3)*(real(Y(6,3)))-E(6)*E(4)*(real(Y(6,4)))-E(6)*E(5)*(real(Y(6,5)))-E(6)*E(7)*(real(Y(6,7)))-E(6)*E(8)*(real(Y(6,8)))-E(6)*E(9)*(real(Y(6,9)))-E(6)*E(10)*(real(Y(6,10))));
        delta_t*(1/H(7,7))*(P(7)-real(Y(7,7))*E(7)^2-E(7)*E(1)*(real(Y(7,1)))-E(7)*E(2)*(real(Y(7,2)))-E(7)*E(3)*(real(Y(7,3)))-E(7)*E(4)*(real(Y(7,4)))-E(7)*E(5)*(real(Y(7,5)))-E(7)*E(6)*(real(Y(7,6)))-E(7)*E(8)*(real(Y(7,8)))-E(7)*E(9)*(real(Y(7,9)))-E(7)*E(10)*(real(Y(7,10))));
        delta_t*(1/H(8,8))*(P(8)-real(Y(8,8))*E(8)^2-E(8)*E(1)*(real(Y(8,1)))-E(8)*E(2)*(real(Y(8,2)))-E(8)*E(3)*(real(Y(8,3)))-E(8)*E(4)*(real(Y(8,4)))-E(8)*E(5)*(real(Y(8,5)))-E(8)*E(6)*(real(Y(8,6)))-E(8)*E(7)*(real(Y(8,7)))-E(8)*E(9)*(real(Y(8,9)))-E(8)*E(10)*(real(Y(8,10))));
        delta_t*(1/H(9,9))*(P(9)-real(Y(9,9))*E(9)^2-E(9)*E(1)*(real(Y(9,1)))-E(9)*E(2)*(real(Y(9,2)))-E(9)*E(3)*(real(Y(9,3)))-E(9)*E(4)*(real(Y(9,4)))-E(9)*E(5)*(real(Y(9,5)))-E(9)*E(6)*(real(Y(9,6)))-E(9)*E(7)*(real(Y(9,7)))-E(9)*E(8)*(real(Y(9,8)))-E(9)*E(10)*(real(Y(9,10))));
        delta_t*(1/H(10,10))*(P(10)-real(Y(10,10))*E(10)^2-E(10)*E(1)*(real(Y(10,1)))-E(10)*E(2)*(real(Y(10,2)))-E(10)*E(3)*(real(Y(10,3)))-E(10)*E(4)*(real(Y(10,4)))-E(10)*E(5)*(real(Y(10,5)))-E(10)*E(6)*(real(Y(10,6)))-E(10)*E(7)*(real(Y(10,7)))-E(10)*E(8)*(real(Y(10,8)))-E(10)*E(9)*(real(Y(10,9))))];
    
    if t<2
        Y = Y3;
    elseif t>=T_fault_in && t<T_fault_end           
        Y = Y32;
    else
        Y = Y3;
    end

    
    yp(1) = delta_t*y(10)+y(1)+u(1);
    yp(2) = delta_t*y(11)+y(2)+u(2);
    yp(3) = delta_t*y(12)+y(3)+u(3);
    yp(4) = delta_t*y(13)+y(4)+u(4);
    yp(5) = delta_t*y(14)+y(5)+u(5);
    yp(6) = delta_t*y(15)+y(6)+u(6);
    yp(7) = delta_t*y(16)+y(7)+u(7);
    yp(8) = delta_t*y(17)+y(8)+u(8);
    yp(9) = delta_t*y(18)+y(9)+u(9);
    yp(10) = delta_t*(1/H(2,2))*(-D(2,2)*(y(10)+u(1))+P(2)-real(Y(2,2))*E(2)^2-...
        E(2)*E(1)*(real(Y(2,1))*cos(y(1))+imag(Y(2,1))*sin(y(1)))-...
        E(2)*E(3)*(real(Y(2,3))*cos(y(1)-y(2))+imag(Y(2,3))*sin(y(1)-y(2)))-...
        E(2)*E(4)*(real(Y(2,4))*cos(y(1)-y(3))+imag(Y(2,4))*sin(y(1)-y(3)))-...
        E(2)*E(5)*(real(Y(2,5))*cos(y(1)-y(4))+imag(Y(2,5))*sin(y(1)-y(4)))-...
        E(2)*E(6)*(real(Y(2,6))*cos(y(1)-y(5))+imag(Y(2,6))*sin(y(1)-y(5)))-...
        E(2)*E(7)*(real(Y(2,7))*cos(y(1)-y(6))+imag(Y(2,7))*sin(y(1)-y(6)))-...
        E(2)*E(8)*(real(Y(2,8))*cos(y(1)-y(7))+imag(Y(2,8))*sin(y(1)-y(7)))-...
        E(2)*E(9)*(real(Y(2,9))*cos(y(1)-y(8))+imag(Y(2,9))*sin(y(1)-y(8)))-...
        E(2)*E(10)*(real(Y(2,10))*cos(y(1)-y(9))+imag(Y(2,10))*sin(y(1)-y(9))))+y(10)-k(1);
    yp(11) = delta_t*(1/H(3,3))*(-D(3,3)*(y(11)+u(2))+P(3)-real(Y(3,3))*E(3)^2-...
        E(3)*E(1)*(real(Y(3,1))*cos(y(2))+imag(Y(3,1))*sin(y(2)))-...
        E(3)*E(2)*(real(Y(3,2))*cos(y(2)-y(1))+imag(Y(3,2))*sin(y(2)-y(1)))-...
        E(3)*E(4)*(real(Y(3,4))*cos(y(2)-y(3))+imag(Y(3,4))*sin(y(2)-y(3)))-...
        E(3)*E(5)*(real(Y(3,5))*cos(y(2)-y(4))+imag(Y(3,5))*sin(y(2)-y(4)))-...
        E(3)*E(6)*(real(Y(3,6))*cos(y(2)-y(5))+imag(Y(3,6))*sin(y(2)-y(5)))-...
        E(3)*E(7)*(real(Y(3,7))*cos(y(2)-y(6))+imag(Y(3,7))*sin(y(2)-y(6)))-...
        E(3)*E(8)*(real(Y(3,8))*cos(y(2)-y(7))+imag(Y(3,8))*sin(y(2)-y(7)))-...
        E(3)*E(9)*(real(Y(3,9))*cos(y(2)-y(8))+imag(Y(3,9))*sin(y(2)-y(8)))-...
        E(3)*E(10)*(real(Y(3,10))*cos(y(2)-y(9))+imag(Y(3,10))*sin(y(2)-y(9))))+y(11)-k(2);
    yp(12) = delta_t*(1/H(4,4))*(-D(4,4)*(y(12)+u(3))+P(4)-real(Y(4,4))*E(4)^2-...
        E(4)*E(1)*(real(Y(4,1))*cos(y(3))+imag(Y(4,1))*sin(y(3)))-...
        E(4)*E(2)*(real(Y(4,2))*cos(y(3)-y(1))+imag(Y(4,2))*sin(y(3)-y(1)))-...
        E(4)*E(3)*(real(Y(4,3))*cos(y(3)-y(2))+imag(Y(4,3))*sin(y(3)-y(2)))-...
        E(4)*E(5)*(real(Y(4,5))*cos(y(3)-y(4))+imag(Y(4,5))*sin(y(3)-y(4)))-...
        E(4)*E(6)*(real(Y(4,6))*cos(y(3)-y(5))+imag(Y(4,6))*sin(y(3)-y(5)))-...
        E(4)*E(7)*(real(Y(4,7))*cos(y(3)-y(6))+imag(Y(4,7))*sin(y(3)-y(6)))-...
        E(4)*E(8)*(real(Y(4,8))*cos(y(3)-y(7))+imag(Y(4,8))*sin(y(3)-y(7)))-...
        E(4)*E(9)*(real(Y(4,9))*cos(y(3)-y(8))+imag(Y(4,9))*sin(y(3)-y(8)))-...
        E(4)*E(10)*(real(Y(4,10))*cos(y(3)-y(9))+imag(Y(4,10))*sin(y(3)-y(9))))+y(12)-k(3);
    yp(13) = delta_t*(1/H(5,5))*(-D(5,5)*(y(13)+u(4))+P(5)-real(Y(5,5))*E(5)^2-...
        E(5)*E(1)*(real(Y(5,1))*cos(y(4))+imag(Y(5,1))*sin(y(4)))-...
        E(5)*E(2)*(real(Y(5,2))*cos(y(4)-y(1))+imag(Y(5,2))*sin(y(4)-y(1)))-...
        E(5)*E(3)*(real(Y(5,3))*cos(y(4)-y(2))+imag(Y(5,3))*sin(y(4)-y(2)))-...
        E(5)*E(4)*(real(Y(5,4))*cos(y(4)-y(3))+imag(Y(5,4))*sin(y(4)-y(3)))-...
        E(5)*E(6)*(real(Y(5,6))*cos(y(4)-y(5))+imag(Y(5,6))*sin(y(4)-y(5)))-...
        E(5)*E(7)*(real(Y(5,7))*cos(y(4)-y(6))+imag(Y(5,7))*sin(y(4)-y(6)))-...
        E(5)*E(8)*(real(Y(5,8))*cos(y(4)-y(7))+imag(Y(5,8))*sin(y(4)-y(7)))-...
        E(5)*E(9)*(real(Y(5,9))*cos(y(4)-y(8))+imag(Y(5,9))*sin(y(4)-y(8)))-...
        E(5)*E(10)*(real(Y(5,10))*cos(y(4)-y(9))+imag(Y(5,10))*sin(y(4)-y(9))))+y(13)-k(4);
    yp(14) = delta_t*(1/H(6,6))*(-D(6,6)*(y(14)+u(5))+P(6)-real(Y(6,6))*E(6)^2-...
        E(6)*E(1)*(real(Y(6,1))*cos(y(5))+imag(Y(6,1))*sin(y(5)))-...
        E(6)*E(2)*(real(Y(6,2))*cos(y(5)-y(1))+imag(Y(6,2))*sin(y(5)-y(1)))-...
        E(6)*E(3)*(real(Y(6,3))*cos(y(5)-y(2))+imag(Y(6,3))*sin(y(5)-y(2)))-...
        E(6)*E(4)*(real(Y(6,4))*cos(y(5)-y(3))+imag(Y(6,4))*sin(y(5)-y(3)))-...
        E(6)*E(5)*(real(Y(6,5))*cos(y(5)-y(4))+imag(Y(6,5))*sin(y(5)-y(4)))-...
        E(6)*E(7)*(real(Y(6,7))*cos(y(5)-y(6))+imag(Y(6,7))*sin(y(5)-y(6)))-...
        E(6)*E(8)*(real(Y(6,8))*cos(y(5)-y(7))+imag(Y(6,8))*sin(y(5)-y(7)))-...
        E(6)*E(9)*(real(Y(6,9))*cos(y(5)-y(8))+imag(Y(6,9))*sin(y(5)-y(8)))-...
        E(6)*E(10)*(real(Y(6,10))*cos(y(5)-y(9))+imag(Y(6,10))*sin(y(5)-y(9))))+y(14)-k(5);
    yp(15) = delta_t*(1/H(7,7))*(-D(7,7)*(y(15)+u(6))+P(7)-real(Y(7,7))*E(7)^2-...
        E(7)*E(1)*(real(Y(7,1))*cos(y(6))+imag(Y(7,1))*sin(y(6)))-...
        E(7)*E(2)*(real(Y(7,2))*cos(y(6)-y(1))+imag(Y(7,2))*sin(y(6)-y(1)))-...
        E(7)*E(3)*(real(Y(7,3))*cos(y(6)-y(2))+imag(Y(7,3))*sin(y(6)-y(2)))-...
        E(7)*E(4)*(real(Y(7,4))*cos(y(6)-y(3))+imag(Y(7,4))*sin(y(6)-y(3)))-...
        E(7)*E(5)*(real(Y(7,5))*cos(y(6)-y(4))+imag(Y(7,5))*sin(y(6)-y(4)))-...
        E(7)*E(6)*(real(Y(7,6))*cos(y(6)-y(5))+imag(Y(7,6))*sin(y(6)-y(5)))-...
        E(7)*E(8)*(real(Y(7,8))*cos(y(6)-y(7))+imag(Y(7,8))*sin(y(6)-y(7)))-...
        E(7)*E(9)*(real(Y(7,9))*cos(y(6)-y(8))+imag(Y(7,9))*sin(y(6)-y(8)))-...
        E(7)*E(10)*(real(Y(7,10))*cos(y(6)-y(9))+imag(Y(7,10))*sin(y(6)-y(9))))+y(15)-k(6);
    yp(16) = delta_t*(1/H(8,8))*(-D(8,8)*(y(16)+u(7))+P(8)-real(Y(8,8))*E(8)^2-...
        E(8)*E(1)*(real(Y(8,1))*cos(y(7))+imag(Y(8,1))*sin(y(7)))-...
        E(8)*E(2)*(real(Y(8,2))*cos(y(7)-y(1))+imag(Y(8,2))*sin(y(7)-y(1)))-...
        E(8)*E(3)*(real(Y(8,3))*cos(y(7)-y(2))+imag(Y(8,3))*sin(y(7)-y(2)))-...
        E(8)*E(4)*(real(Y(8,4))*cos(y(7)-y(3))+imag(Y(8,4))*sin(y(7)-y(3)))-...
        E(8)*E(5)*(real(Y(8,5))*cos(y(7)-y(4))+imag(Y(8,5))*sin(y(7)-y(4)))-...
        E(8)*E(6)*(real(Y(8,6))*cos(y(7)-y(5))+imag(Y(8,6))*sin(y(7)-y(5)))-...
        E(8)*E(7)*(real(Y(8,7))*cos(y(7)-y(6))+imag(Y(8,7))*sin(y(7)-y(6)))-...
        E(8)*E(9)*(real(Y(8,9))*cos(y(7)-y(8))+imag(Y(8,9))*sin(y(7)-y(8)))-...
        E(8)*E(10)*(real(Y(8,10))*cos(y(7)-y(9))+imag(Y(8,10))*sin(y(7)-y(9))))+y(16)-k(7);
    yp(17) = delta_t*(1/H(9,9))*(-D(9,9)*(y(17)+u(8))+P(9)-real(Y(9,9))*E(9)^2-...
        E(9)*E(1)*(real(Y(9,1))*cos(y(8))+imag(Y(9,1))*sin(y(8)))-...
        E(9)*E(2)*(real(Y(9,2))*cos(y(8)-y(1))+imag(Y(9,2))*sin(y(8)-y(1)))-...
        E(9)*E(3)*(real(Y(9,3))*cos(y(8)-y(2))+imag(Y(9,3))*sin(y(8)-y(2)))-...
        E(9)*E(4)*(real(Y(9,4))*cos(y(8)-y(3))+imag(Y(9,4))*sin(y(8)-y(3)))-...
        E(9)*E(5)*(real(Y(9,5))*cos(y(8)-y(4))+imag(Y(9,5))*sin(y(8)-y(4)))-...
        E(9)*E(6)*(real(Y(9,6))*cos(y(8)-y(5))+imag(Y(9,6))*sin(y(8)-y(5)))-...
        E(9)*E(7)*(real(Y(9,7))*cos(y(8)-y(6))+imag(Y(9,7))*sin(y(8)-y(6)))-...
        E(9)*E(8)*(real(Y(9,8))*cos(y(8)-y(7))+imag(Y(9,8))*sin(y(8)-y(7)))-...
        E(9)*E(10)*(real(Y(9,10))*cos(y(8)-y(9))+imag(Y(9,10))*sin(y(8)-y(9))))+y(17)-k(8);
    yp(18) = delta_t*(1/H(10,10))*(-D(10,10)*(y(18)+u(9))+P(10)-real(Y(10,10))*E(10)^2-...
        E(10)*E(1)*(real(Y(10,1))*cos(y(9))+imag(Y(10,1))*sin(y(9)))-...
        E(10)*E(2)*(real(Y(10,2))*cos(y(9)-y(1))+imag(Y(10,2))*sin(y(9)-y(1)))-...
        E(10)*E(3)*(real(Y(10,3))*cos(y(9)-y(2))+imag(Y(10,3))*sin(y(9)-y(2)))-...
        E(10)*E(4)*(real(Y(10,4))*cos(y(9)-y(3))+imag(Y(10,4))*sin(y(9)-y(3)))-...
        E(10)*E(5)*(real(Y(10,5))*cos(y(9)-y(4))+imag(Y(10,5))*sin(y(9)-y(4)))-...
        E(10)*E(6)*(real(Y(10,6))*cos(y(9)-y(5))+imag(Y(10,6))*sin(y(9)-y(5)))-...
        E(10)*E(7)*(real(Y(10,7))*cos(y(9)-y(6))+imag(Y(10,7))*sin(y(9)-y(6)))-...
        E(10)*E(8)*(real(Y(10,8))*cos(y(9)-y(7))+imag(Y(10,8))*sin(y(9)-y(7)))-...
        E(10)*E(9)*(real(Y(10,9))*cos(y(9)-y(8))+imag(Y(10,9))*sin(y(9)-y(8))))+y(18)-k(9);

end

%------------- END OF CODE --------------