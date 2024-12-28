% Demo 3
% Compute various functions 
% of the asymptotic covariance matrix
% in particular 
% a/ det(V)
% b/ trace(V)
% c/ condition_number(V)
% d/ max(eigenvalue(V))

clear all
close all
clc


% define class
% use sensitivity class
a = ecr_sens;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Class Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

% input into class the value of k and D
a.k_ref = 1.5E-4; % unit of k is m/s
a.D_ref = 2.5E-9; % unit of D m^2/s

a.standard_dev = 1E-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Sigma
% and the Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Timespan
t_min = 0.0;
t_max = 1E4;

time = linspace(t_min, t_max, 1E4+1)'; % Maximum time is 10,000 s

% vector of half-sizes
half_thickness_min = 1E-4; % 1E-4 m = 100 \mu m
half_thickness_max = 1E-2; % 1E-2 m = 1 cm
half_thickness_vec = logspace(log10(half_thickness_min), log10(half_thickness_max), 100);

det_cov_sigma_n_vec = zeros(size(half_thickness_vec));
trace_cov_sigma_n_vec = zeros(size(half_thickness_vec));
cond_cov_sigma_n_vec = zeros(size(half_thickness_vec));
max_eig_cov_sigma_n_vec = zeros(size(half_thickness_vec));

for iter_half_thickness = 1 : numel(half_thickness_vec)

    % input into class the size of the sample
    a.half_thickness_x_ref = half_thickness_vec(iter_half_thickness); % unit of x-thickness is m
    a.half_thickness_y_ref = half_thickness_vec(iter_half_thickness); % unit of y-thickness is m
    a.half_thickness_z_ref = half_thickness_vec(iter_half_thickness); % unit of z-thickness is m
    
    fprintf('iteration = %u/%u \n', iter_half_thickness, numel(half_thickness_vec))
    det_cov_sigma_n_vec(iter_half_thickness) = a.det_cov_sigma_n(time);       
    trace_cov_sigma_n_vec(iter_half_thickness) = a.trace_cov_sigma_n(time);
    cond_cov_sigma_n_vec(iter_half_thickness) = a.cond_cov_sigma_n(time);
    max_eig_cov_sigma_n_vec(iter_half_thickness) = a.max_eig_cov_sigma_n(time);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(2,2,1)
loglog(half_thickness_vec, det_cov_sigma_n_vec, '-k', 'LineWidth', 3)
xlabel('thickness/m', 'Interpreter', 'Latex', 'FontSize', 30)
ylabel('det(V)', 'Interpreter', 'Latex', 'FontSize', 30)
set(gca,'FontSize',20)

subplot(2,2,2)
loglog(half_thickness_vec, trace_cov_sigma_n_vec, '-k', 'LineWidth', 3)
xlabel('thickness/m', 'Interpreter', 'Latex', 'FontSize', 30)
ylabel('tr(V)', 'Interpreter', 'Latex', 'FontSize', 30)
set(gca,'FontSize',20)

subplot(2,2,3)
loglog(half_thickness_vec, cond_cov_sigma_n_vec, '-k', 'LineWidth', 3)
xlabel('thickness/m', 'Interpreter', 'Latex', 'FontSize', 30)
ylabel('$\kappa(\rm V)$', 'Interpreter', 'Latex', 'FontSize', 30)
set(gca,'FontSize',20)

subplot(2,2,4)
loglog(half_thickness_vec, max_eig_cov_sigma_n_vec, '-k', 'LineWidth', 3)
xlabel('thickness/m', 'Interpreter', 'Latex', 'FontSize', 30)
ylabel('$\max(\lambda(\rm V))$', 'Interpreter', 'Latex', 'FontSize', 30)
set(gca,'FontSize',20)

X = 45.0;                  % paper size
Y = 30.0;                  % A3 paper size
xMargin = 3.;              % left/right margins from page borders
yMargin = 3.;              % bottom/top margins from page borders
xSize = X + 2*xMargin;     %# figure size on paper (width & height)
ySize = Y + 2*yMargin;     %# figure size on paper (width & height)

%# figure size on screen (50% scaled, but same aspect ratio)
set(gcf, 'Units','centimeters', 'Position',[3 3 xSize ySize]/2)
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')
