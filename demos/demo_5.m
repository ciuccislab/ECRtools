% demo_5
% In this demo we will:
% 
% 1. perform synthetic experiments 
% 2. fit the synthetic experiments
% 3. plot the outcome of the synthetic experiment
%       fits and output that to figure along with the
%       asymptotic covariance matrix
%
% the steps are described in detail below


clear all
close all
clc

% define experimental time
t_exp = linspace(0, 500, 151)';

% declare the fitting class, which has ecr_sens (ecr sensitivity class)
% as superclass and upload the parameters into the class
a = ecr_fit;
a.k_ref = 1E-4; % k - unit of k is m/s
a.D_ref = 1E-9; % D - unit of D m^2/s
a.half_thickness_x_ref = 1E-3; % half thickness - unit is m
a.half_thickness_y_ref = 1E-3; % half thickness - unit is m
a.half_thickness_z_ref = 1E-3; % half thickness - unit is m
a.standard_dev = 1E-2; % 1% error
a.meas_time = t_exp; % experimental time - unit is s

% note that the measurement itself is not uploaded as it will change
% through the computation

% ouput the boundary of the covariance matrix
% (theta)^T*V^-1*(theta)<3^2 (this number)
boundary_cov = a.boundary_cov_sigma_n(a.meas_time);

% run synthetic experiment
tic 
% note that one can use
% several optimizers by 
% theta_vec = a.synthetic_exp('pswarm', 1000);
% theta_vec = a.synthetic_exp('nomad', 1000);
% theta_vec = a.synthetic_exp('GN_DIRECT', 1000);
theta_vec = a.synthetic_exp('nlopt', 100);
% note that by selecting different internal optimizer
% the time of the class call changes
toc

figure(1)
% plot of the covariance matrix
plot(boundary_cov(1,:), boundary_cov(2,:), '-k', 'LineWidth', 2)
hold on
xlabel('$\hat k/k_{\rm exact}$', 'Interpreter', 'Latex', 'FontSize', 25)
ylabel('$\hat D/D_{\rm exact}$', 'Interpreter', 'Latex', 'FontSize', 25)
set(gca,'FontSize',14)
axis([.4, 1.6, 0.4, 1.6]);
% plot of the covariance matrix
plot(theta_vec(1,:), theta_vec(2,:), '.r', 'LineWidth', 2)

X = 35.0;                  % paper size
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

clearvars -except a