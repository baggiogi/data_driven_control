function [control_energy_mb, control_energy_dd, err_mb, err_dd] = data_driven_fMRI(subject_id)
% DATA_DRIVEN_FMRI - Computes the control energies and errors in the final
% state for the data-driven input and model-based input
%
% Task-fMRI has been data extracted from the HCP S1200 dataset:
%
% https://www.humanconnectome.org/study/hcp-young-adult/document/1200-subjects-data-release
%
% The model-based input is computed using a linear model estimated as in:
%
% C. O. Becker, D. S. Bassett, V. M. Preciado.
% "Large-scale dynamic modeling of task-fMRI signals via subspace system identification."
% Journal of neural engineering 15.6 (2018): 066016.
%
% Syntax:  [control_energy_mb, control_energy_dd, err_mb, err_dd] = data_driven_fMRI(subject_id)
%
% Inputs:
%    subject_id - subject id of HCP dataset
%
% Outputs:
%    control_energy_mb - control energy model-based input
%    control_energy_dd - control energy data-driven input
%    err_mb - error of model-based input
%    err_dd - error of data-driven input
%
% Author: Giacomo Baggio
% email: baggio@dei.unipd.it
% March 2020; Last revision: 01-April-2020

%------------- BEGIN CODE --------------

%% Parameters
r = 148;        % number of brain regions
L = 284;        % total time horizon
T = 100;        % control horizon
T_est = 150;    % estimation horizon
m = 6;          % number of inputs
p = 148;        % number of outputs

% filter design parameters
A_stop1 = 20;		% Attenuation in the first stopband [Hz]
F_stop1 = 0.04;		% Edge of the stopband [Hz]
F_pass1 = 0.06; 	% Edge of the passband [Hz]
F_pass2 = 0.12;     % Closing edge of the passband [Hz]
F_stop2 = 0.15;     % Edge of the second stopband [Hz]
A_stop2 = 20;		% Attenuation in the second stopband [Hz]
A_pass  = 1;        % Amount of ripple allowed in the passband [Hz]

%% Construct data matrices

% data for data-driven control
U   = [];
Yf  = [];

% data for system identification
Ui = [];
Yi = [];

% load data of subject subject_id
BOLD    = ['./fMRI_data/bold_' num2str(subject_id) '.txt'];
CUES    = ['./fMRI_data/cues_' num2str(subject_id) '.txt'];
HEART   = ['./fMRI_data/heart_' num2str(subject_id) '.txt'];
RESP    = ['./fMRI_data/resp_' num2str(subject_id) '.txt'];

% crop data to interval [1,284]
outputs = load(BOLD);
outputs = outputs(:,1:284);
inputs  = load(CUES);
heart   = load(HEART)';
heart   = heart(1:284);
resp    = load(RESP)';
resp    = resp(1:284);

% encode visual cues
inputs(6,12)  = 1; inputs(6,33)  = 1; inputs(6,54)  = 1; inputs(6,75)  = 1;
inputs(6,96)  = 1; inputs(6,138) = 1; inputs(6,159) = 1; inputs(6,180) = 1;
inputs(6,221) = 1; inputs(6,242) = 1;

% rearrange inputs: 1) CUE 2) LF 3) LH 4) RF 5) RH 6) T
inputs = inputs([6 4 3 2 5 1],1:284);

% physiological signals
U_P = [heart; resp];

% physiologically regressed outputs
outputs_r = outputs*(eye(L)-pinv(U_P)*U_P);

% band-pass filter
Hd = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',50,F_stop1,F_pass1,F_pass2,F_stop2);
d  = design(Hd,'equiripple');

% filtered outputs
outputs_rf = filter(d,outputs_r');
outputs_rf = outputs_rf';

inputs_tmp = [zeros(m,T-10),inputs];
outputs_tmp = [zeros(p,T-10),outputs_rf];

Yi = reshape(outputs_rf(:,1:T_est),[],1);
Ui = reshape(inputs(:,1:T_est),[],1);

for t_s =1:T
    inputs_sel = flip(inputs_tmp(:,t_s:t_s+T-1),2);
    inputs_r = reshape(inputs_sel,[],1);
    U = [U inputs_r];
    Yf = [Yf outputs_tmp(:,t_s+T)];
end

%% Model identification

% subspace identification parameters
s = 3;
len = T_est+1-s;
Yh = zeros(s*p,len);
Uh = zeros(s*m,len);

% convert data in Hankel form
Y_tmp = Yi;
U_tmp = Ui;

for j = 1:len
    Yh(:,j) = Y_tmp(1+p*(j-1):p*(j-1)+s*p);
    Uh(:,j) = U_tmp(1+m*(j-1):m*(j-1)+s*m);
end

Uperp = eye(len) - Uh'*pinv(Uh*Uh')*Uh;
M = Yh*Uperp;
[Us,Sigma,Vs] = svd(M);
% select model order
nr = 20;
%nr = max([20,min(find(diag(Sigma/len)<0.0115))])
% extended observability estimate
Un = Us(:,1:nr);

% estimate A matrix
Aest = pinv(Un(1:p*(s-1),1:nr))*Un((p+1):p*s,1:nr);

% enforce stability of A
if max(abs(eig(Aest)))>1
    Aest = Aest/(max(abs(eig(Aest)))+1e-2);
end

% estimate C matrix
Cest = Un(1:p,1:nr);

% compute matrix K
  
Y_tmp = Yi;
U_tmp = Ui;
K = [];

for j = 1:T_est
    
    K_tmp = zeros(p,nr*m);
    
    if j ~= 1
        for r = 2:j
            K_tmp = K_tmp+kron(U_tmp(m*(r-2)+1:m*(r-1))',Cest*Aest^(j-r));
        end
    end
    
    K = [K; K_tmp];
end

Yh_tot = Y_tmp;

% estimate B matrix
gamma = 5; % regularization parameter
Best = inv(K'*K+gamma*eye(m*nr))*K'*Yh_tot;
Best = reshape(Best,[nr,m]);

% estimated system
sys_est = ss(Aest,Best,Cest,[],-1);


%% compute error and control energy

% compute output controllability matrix
Co = sys_est.C*sys_est.B;

for jj=1:T-1
    Co = [Co sys_est.C*sys_est.A^jj*sys_est.B];
end

% set of target states
reach_outputs = orth(Yf*pinv(U));
n_reach = nr;

control_energy_mb = zeros(n_reach,1);
control_energy_dd = zeros(n_reach,1);
err_mb = zeros(n_reach,1);
err_dd = zeros(n_reach,1);

for j=1:n_reach
    
    yf = reach_outputs(:,j);
    
    % model-based control
    u_mb = Co'*pinv(Co*Co')*yf;
    % data-driven control
    u_dd = pinv(Yf*pinv(U))*yf;
    
    y_mb = Co*u_mb;
    y_dd = Co*u_dd;
    
    % normalized error model-based
    err_mb(j) = norm(yf-y_mb)/p;
    % normalized error data-driven
    err_dd(j) = norm(yf-y_dd)/p;
    
    % control energy model-based
    control_energy_mb(j) = norm(u_mb);
    % control energy data-driven
    control_energy_dd(j) = norm(u_dd);
    
end

%------------- END OF CODE --------------