% demo_4
% In this demo we will:
%
% 1. fit experimental data
% 2. plot the outcome of the fit by showing simultaneously
%       experimental data, fit and residual
% 3. plot the asymptotic confidence region
% the steps are described in detail below

clear all
close all
clc

k_stored = [];
D_stored = [];

for iter = 1: 20
    
    ECR_data = load('data.txt');
    time_exp = ECR_data(:,1)-ECR_data(1,1); % measurement_time
    sigma_n_exp = ECR_data(:,2); % measured \sigma_n
    
    a = ecr_fit;
    
    % input into class the value of k and D
    
    
    a.k_ref = abs(1E-4*(1+0.2*randn)); % unit of k is m/s
    a.D_ref = abs(1E-6*(1+0.2*randn)); % unit of D m^2/s
    
    % input into class the half size of the sample
    a.half_thickness_x_ref = 1E-3; % unit of half x-thickness is m
    a.half_thickness_y_ref = 1E-3; % unit of half y-thickness is m
    a.half_thickness_z_ref = 1E-3; % unit of half z-thickness is m
    
    % upload the error of made in the measurement:
    % standard deviation = 0.01
    a.standard_dev = 1E-2; % this will be recomputed later
    
    a.meas_time = time_exp;
    a.meas_sigma_n = sigma_n_exp;
    
    % first run fit within a large span
    theta_out = a.log10_fit();
    
    % update the reference value
    a.k_ref = theta_out(1)*a.k_ref;
    a.D_ref = theta_out(2)*a.D_ref;
    
    theta_out = a.fit('nomad');
    a.k_ref = theta_out(1)*a.k_ref;
    a.D_ref = theta_out(2)*a.D_ref;
    
    sigma_computed = a.output_sigma_n_det(theta_out);
    error_meas = (sigma_computed-sigma_n_exp);
    std_computed = a.compute_std_meas();
    a.standard_dev = std_computed;
    
    k_stored = [k_stored; a.k_ref];
    D_stored = [D_stored; a.D_ref];
    
    figure(1)
    
    [AX,H1,H2] = plotyy(time_exp, sigma_computed, time_exp, 100*error_meas,'plot','plot');
    hold on
    
    set(H1,'LineStyle','-', 'Linewidth', 3, 'Color', 'k');
    set(H2,'Marker','+', 'MarkerSize', 4.5, 'Color', 'b');
    
    set(AX(1),'ycolor','k')
    set(AX(1),'FontSize',15)
    set(AX(1),'XLim',[0 max(time_exp)])
    set(AX(1),'YLim',[0 1.0])
    set(AX(1),'YTick',[0:0.1:1.0])
    set(AX(2),'ycolor','b')
    set(AX(2),'FontSize',15)
    set(AX(2),'XLim',[0 max(time_exp)])
    set(AX(2),'YLim',[-100 100])
    set(AX(2),'YTick',[-100:10:100])
    
    plot(time_exp, sigma_n_exp, '+r', 'MarkerSize', 4.5)
    plot([0, max(time_exp)], [3*std_computed, 3*std_computed]+0.5, '--g', 'Linewidth', 3);
    plot([0, max(time_exp)], [-3*std_computed, -3*std_computed]+0.5, '--g', 'Linewidth', 3);
    
    hold off
    xlabel('$t_{\rm exp}/{\rm s}$', 'Interpreter' , 'Latex', 'Fontsize', 25)
    %set(get(AX(1),'Ylabel'),'String','$\sigma_n$', 'Interpreter' , 'Latex', 'Fontsize', 18)
    set(get(AX(1),'Ylabel'),'String','$\sigma_n$', 'Interpreter', 'Latex', 'Fontsize', 25)
    set(get(AX(2),'Ylabel'),'String','$100\times(\sigma_n-\sigma_n^{\rm meas}~)$', 'Interpreter' , 'Latex', 'Fontsize', 25)
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
    
    clearvars -except a k_stored D_stored
    
    % a further capability of ECR is that
    % it can output confidence ellipsoids
    boundary_cov = a.boundary_cov_sigma_n(a.meas_time);
    
    figure(2)
    plot(boundary_cov(1,:), boundary_cov(2,:), '-k', 'LineWidth', 2)
    xlabel('$\hat k/k_{\rm fit}$', 'Interpreter', 'Latex', 'FontSize', 25)
    ylabel('$\hat D/D_{\rm fit}$', 'Interpreter', 'Latex', 'FontSize', 25)
    set(gca,'FontSize',14)
    axis([0, 10, 0.0, 10]);
    
    % plot of the covariance matrix
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
    
    clearvars -except k_stored D_stored

end

figure(3)
hist(log10(D_stored));
xlabel('log10(D)')

figure(4)
hist(log10(k_stored));
xlabel('log10(k)')