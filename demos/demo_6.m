% demo_6
% In this demo we will optimize experiments
% provided that we have a value for 
% the measurement error
% and for the experimental parameters k and D
%
% the steps described in detail below

clear all
close all
clc

% declare the fitting class, which has ecr_sens (ecr sensitivity class)
% as superclass and upload the parameters into the class
a = ecr_opti;
a.k_ref = 1E-4; % k - unit of k is m/s
a.D_ref = 1E-9; % D - unit of D m^2/s
a.half_thickness_x_ref = 1E-3; % half thickness - unit is m
a.half_thickness_y_ref = 1E-3; % half thickness - unit is m
a.half_thickness_z_ref = 1E-3; % half thickness - unit is m
a.standard_dev = 1E-2; % 1% error

% define the bounds for the experimental design
% what time span is acceptable?
time_exp_bound = [1E2, 1E4];
% which thicknesses are ok?
half_thickness_bound = [5E-4, 5E-2];
% fix the number of measurements per second
N_meas_per_second = 1;
% establish the criterion used for optimization
% in this case, select the determinant
criterion = 'det';
% we define the criterion and optimize 
% _Ntot = number of measurement points per second
[time_det_opt, thickness_det_opt, det_opt] = a.optimize_exp_Nps(time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, '3D');

% I show below other options
% criterion = 'trace';
% [time_trace_opt, thickness_trace_opt, trace_opt] = a.optimize_exp_Nps(time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, '3D');
% criterion = 'max_eig';
% [time_max_eig_opt, thickness_max_eig_opt, max_eig_opt] = a.optimize_exp_Nps(time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, '3D');
% criterion = 'cond';
% [time_cond_opt, thickness_cond_opt, cond_opt] = a.optimize_exp_Nps(time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, '3D');

% if instead the total number of measurements 
% is fixed rather than the number of measurements per second
% then we declare N_tot (number of time points)
N_tot = 1E2;
% we define the criterion and optimize _Ntot = number total points
criterion = 'det';
[time_det_opt2, thickness_det_opt2, det_opt2] = a.optimize_exp_Ntot(time_exp_bound, half_thickness_bound, N_tot, criterion, '3D');
% criterion = 'trace';
% [time_trace_opt2, thickness_trace_opt2, trace_opt2] = a.optimize_exp_Ntot(time_exp_bound, half_thickness_bound, N_tot, criterion, '3D');
% criterion = 'max_eig';
% [time_max_eig_opt2, thickness_max_eig_opt2, max_eig_opt2] = a.optimize_exp_Ntot(time_exp_bound, half_thickness_bound, N_tot, criterion, '3D');
% criterion = 'cond';
% [time_cond_opt2, thickness_cond_opt2, cond_opt2] = a.optimize_exp_Ntot(time_exp_bound, half_thickness_bound, N_tot, criterion, '3D');


% below we check if the results obtained are reasonable 
% by computing the two objective functions
% over the admissible design set

time_exp_vec = logspace(log10(min(time_exp_bound)), log10(max(time_exp_bound)), 11);
half_thickness_vec = logspace(log10(min(half_thickness_bound)), log10(max(half_thickness_bound)), 11);

[time_exp_mat, half_thickness_mat] = meshgrid(time_exp_vec, half_thickness_vec);

[number_rows, number_cols] = size(time_exp_mat);

det_cov_sigma_n = zeros(size(time_exp_mat));
det_cov_sigma_n2 = det_cov_sigma_n;

for iter_rows = 1 : number_rows
    
    for iter_cols = 1 : number_cols
        
        a.half_thickness_x_ref = half_thickness_mat(iter_rows, iter_cols); % half thickness - unit is m
        a.half_thickness_y_ref = half_thickness_mat(iter_rows, iter_cols); % half thickness - unit is m
        a.half_thickness_z_ref = half_thickness_mat(iter_rows, iter_cols); % half thickness - unit is m
        
        % this needs to be a vector!!
        time    = [0: 1/N_meas_per_second: time_exp_mat(iter_rows, iter_cols)]';
        time2   = linspace(0, time_exp_mat(iter_rows, iter_cols),N_tot)';
        
        det_cov_sigma_n(iter_rows, iter_cols) = a.det_cov_sigma_n(time); 
        det_cov_sigma_n2(iter_rows, iter_cols) = a.det_cov_sigma_n(time2); 

        fprintf('(%u, %u), thickness = %e, time_max = %e, det(V) = %e, det(V2) = %e\n',...
            iter_rows, iter_cols, ...
            half_thickness_mat(iter_rows, iter_cols), time_exp_mat(iter_rows, iter_cols), ...
            det_cov_sigma_n(iter_rows, iter_cols), det_cov_sigma_n2(iter_rows, iter_cols));
        
    end
    
end

% then plot the determinant of the covariance matrix

figure(1)
h = surf(time_exp_mat, half_thickness_mat, det_cov_sigma_n)
hold on 
plot3(time_det_opt, thickness_det_opt, det_opt, '.r', 'MarkerSize', 40)
set(get(h,'Parent'),'XScale','log');
set(get(h,'Parent'),'YScale','log');
set(get(h,'Parent'),'ZScale','log');
xlabel('t_{exp}/s')
ylabel('L/m')
zlabel('\phi_D')
view([65,25])
axis([1E2, 1E4, 1E-4, 1E-1, 1E-10, 1E0])

figure(2)
h = surf(time_exp_mat, half_thickness_mat, det_cov_sigma_n2)
hold on 
plot3(time_det_opt2, thickness_det_opt2, det_opt2, '.r', 'MarkerSize', 40)
set(get(h,'Parent'),'XScale','log');
set(get(h,'Parent'),'YScale','log');
set(get(h,'Parent'),'ZScale','log');
xlabel('t_{exp}/s')
ylabel('L/m')
zlabel('\phi_D')
axis([1E2, 1E4, 1E-4, 1E-1, 1E-10, 1E0])
