% PLOT_RESULTS_FMRI - Plots average error and control energy of model-based 
% and data-driven controls of task-based fMRI dynamics
%
% Other m-files required: data_driven_fMRI.m
%
% Author: Giacomo Baggio
% email: baggio@dei.unipd.it
% March 2020; Last revision: 01-April-2020

%------------- BEGIN CODE --------------

clear all
close all
clc

% number of subjects
N = 1;          
% number of target states
nr = 20;

% load task-fMRI data
load ./fMRI_data/unrelated_subjects_final.txt
% pick N subjects randomly 
subjects_tot = datasample(1:280,N,'Replace',false);

err_mb_tot = zeros(nr,N);
err_dd_tot = zeros(nr,N);
control_energy_mb_tot = zeros(nr,N);
control_energy_dd_tot = zeros(nr,N);

%% compute errors and control energies
for ii = 1:N
    
    subject_id = unrelated_subjects_final(subjects_tot(ii));
    [control_energy_mb, control_energy_dd, err_mb, err_dd] = data_driven_fMRI(subject_id);
    
    err_mb_tot(:,ii) = err_mb;
    err_dd_tot(:,ii) = err_dd;
    control_energy_mb_tot(:,ii) = control_energy_mb;
    control_energy_dd_tot(:,ii) = control_energy_dd;
    
end

%% plot results
figure
subplot(2,1,1)

mean_err_mb = sum(err_mb_tot,2)/N;
std_err_mb = std(err_mb_tot,[],2);
mean_err_dd = sum(err_dd_tot,2)/N;
std_err_dd = std(err_dd_tot,[],2);
mean_control_energy_mb = sum(control_energy_mb_tot,2)/N;
std_control_energy_mb = std(control_energy_mb_tot,[],2);
mean_control_energy_dd = sum(control_energy_dd_tot,2)/N;
std_control_energy_dd = std(control_energy_dd_tot,[],2);

bar([mean_err_mb mean_err_dd])
legend('model-based','data-driven')
title('norm of the error on final output')
subplot(2,1,2)
bar([mean_control_energy_mb mean_control_energy_dd])
title('control energy')
legend('model-based','data-driven')

%------------- END OF CODE --------------
