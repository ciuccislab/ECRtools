clear all
close all
clc

% Demo 1
% Computed Normalized Conductivity
% and Synthetic Measurement

% define class
a = ecr_base;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Class Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%

% input into class the value of k and D
a.k_ref = 1.5E-4; % unit of k is m/s
a.D_ref = 2.5E-9; % unit of D m^2/s

% input into class the half size of the sample
a.half_thickness_x_ref = 1E-3; % unit of x-thickness is m
a.half_thickness_y_ref = 1E-3; % unit of y-thickness is m
a.half_thickness_z_ref = 1E-3; % unit of z-thickness is m

% upload the error of made in the measurement: 
% standard deviation = 0.01
a.standard_dev = 1E-2;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Sigma
% and the Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Timespan
N = 101;
time = linspace(0, 1E2, N)'; % Maximum time is 100 s

% compute the exact normalized conductivity \sigma_n
sigma_n = a.sigma_n_det(time);

% synthetic measurement
sigma_n_meas = a.sigma_n_meas(time);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(time, sigma_n, '-k', 'LineWidth', 3)
hold on
plot(time, sigma_n_meas, '.r', 'MarkerSize', 15)

xlabel('$t$/s', 'Interpreter', 'Latex', 'FontSize', 30)
ylabel('$\sigma_n$', 'Interpreter', 'Latex', 'FontSize', 30)
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
