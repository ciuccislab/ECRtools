% demo_7
%
% In this demo we will optimize experiments
% provided that we have an estimate 
% for the measurement error
% and a range for the experimental parameters k and D
%
% the steps described in detail below

clear all
close all
clc

% declare the optimization class, which has ecr_fit (ecr fitting class)
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
time_exp_bound = [5E2, 5E4];

% which thicknesses are ok?
half_thickness_bound = [5E-4, 5E-2];

% expected span for the robust design
k_int = linspace(0.1*a.k_ref,a.k_ref, 11);
D_int = linspace(a.D_ref,10.*a.D_ref, 11);


% establish the criterion used for optimization
% in this case, select the determinant
criterion = 'det';
% we define the criterion and optimize 
% _Ntot = number of measurement points per second
% then we declare N_tot (number of time points)


% we optimize for _Ntot = number total points
N_tot = 151;
[time_opt_avg_Ntot, thickness_opt_avg_Ntot, goal_opt] = a.optimize_avg_robust_Ntot(time_exp_bound, half_thickness_bound, N_tot, criterion, '3D', k_int, D_int);
[time_opt_minmax_Ntot, thickness_opt_minmax_Ntot, goal_opt] = a.optimize_minmax_robust_Ntot(time_exp_bound, half_thickness_bound, N_tot, criterion, '3D', k_int, D_int);

% we optimize for _Nps = number of points per second
% fix the number of measurements per second
N_meas_per_second = 1;
[time_opt_avg_Nps, thickness_opt_avg_Nps, goal_opt] = a.optimize_avg_robust_Nps(time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, k_int, D_int);
[time_opt_minmax_Nps, thickness_opt_minmax_Nps, goal_opt] = a.optimize_minmax_robust_Nps(time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, k_int, D_int);      