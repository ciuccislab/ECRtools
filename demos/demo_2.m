clear all
close all
clc

% Demo 2
% Compute Sensitivitites
% of the Measurement 
% with Respect to the k and D

% define class
% use sensitivity class
a = ecr_sens;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Sigma
% and the Measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define timespan for gradient computation
N = 1000;
t_min = 1E-2; % Minimum time is 1E-2 s
t_max = 1E4; % Maximum time is 1E4 s
time = logspace(log10(t_min), log10(t_max), 1000)';

grad_out = a.grad_sigma_n(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
semilogx(time, grad_out(:, 1), '-k', 'LineWidth', 3)
hold on
semilogx(time, grad_out(:, 2), '-r', 'LineWidth', 3)

xlabel('$t$/s', 'Interpreter', 'Latex', 'FontSize', 30)
ylabel('$q\frac{\partial \sigma_n}{\partial q}$', 'Interpreter', 'Latex', 'FontSize', 30)
set(gca,'FontSize',20)
legend('q=k', 'q=D')

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
