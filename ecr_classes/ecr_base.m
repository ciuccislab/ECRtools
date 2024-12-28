% This file is part of ECR_toolbox - Electrical Conductivity Relaxation Toolbox ver 0.01.
%
% Copyright (C)2012  Francesco Ciucci
%
% Department of Mechanical Engineering 
% The Hong Kong University of Science and Technology
% Clearwater Bay
% Hong Kong, China SAR
%
% Send bug reports and feedback to: mefrank@ust.hk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


classdef ecr_base
    
    % base parameters for ECR_output
    properties
    
        k_ref % reference value for k
        D_ref % reference value for D

 
        half_thickness_x_ref % reference half thickness in x direction
        half_thickness_y_ref % reference half thickness in y direction
        half_thickness_z_ref % reference half thickness in z_direction

        standard_dev % standard deviation
    
    end
   
    % base private parameters for ECR_output
%     properties (SetAccess = private, GetAccess = private)
%         
%         theta_k = 1.0 % theta_k = k/k_ref
%         theta_D = 1.0 % theta_D = D/D_ref
%         
%         xi_x = 1.0 % x/x_ref
%         xi_y = 1.0 % y/y_ref
%         xi_z = 1.0 % z/z_ref
%         
%     end
    
    % methods/functions in the class accessible to the user
    methods

        
        function out_sigma_n = sigma_n_det(obj, time)

            half_thickness_x = obj.half_thickness_x_ref;
            half_thickness_y = obj.half_thickness_y_ref;
            half_thickness_z = obj.half_thickness_z_ref;

            if half_thickness_x == Inf || half_thickness_y == Inf || half_thickness_z == Inf
                if half_thickness_x == Inf && half_thickness_y<Inf && half_thickness_z<Inf
                    out_sigma_n = sigma_n_det_2D(obj, time);
                elseif half_thickness_y == Inf && half_thickness_x<Inf && half_thickness_z<Inf
                    out_sigma_n = sigma_n_det_2D(obj, time);
                elseif half_thickness_z == Inf && half_thickness_x<Inf && half_thickness_y<Inf
                    out_sigma_n = sigma_n_det_2D(obj, time);
                else
                    out_sigma_n = sigma_n_det_1D(obj, time);
                end
            elseif half_thickness_x == half_thickness_y && half_thickness_x == half_thickness_z
                out_sigma_n = sigma_n_det_symm(obj, time);
            else
                out_sigma_n = sigma_n_det_unsymm(obj, time);
            end
        end
        
        function out_sigma_n_meas = sigma_n_meas(obj, time)
            out_sigma_n_det = sigma_n_det(obj, time);
            out_error_sigma = error_fct(obj, out_sigma_n_det);
            out_sigma_n_meas = out_sigma_n_det + out_error_sigma;
        end
        
    end
    
    % methods/functions in the class not accessible to the user
    methods (Access = protected)
        
        % sigma_n for cube
        function out_sigma_n = sigma_n_det_symm(obj, time)

            k = obj.k_ref;
            D = obj.D_ref;

            half_thickness_x = obj.half_thickness_x_ref;
            %half_thickness_y = xi_y*half_thickness_y_ref;
            %half_thickness_z = xi_z*half_thickness_z_ref;

            L_x = half_thickness_x*k/D;
            %L_y = half_thickness_y*k/D;
            %L_z = half_thickness_z*k/D;

            beta_x = obj.diffbeta(L_x);
            %beta_y = beta_x; diffbeta(L_y);
            %beta_z = beta_x; diffbeta(L_z);

            out_sigma_n = zeros(size(time));
            coeff_x = 2*L_x^2./(beta_x.^2.*(beta_x.^2+L_x^2+L_x));
            %coeff_y = 2*L_y^2./(beta_y.^2.*(beta_y.^2+L_y^2+L_y));
            %coeff_z = 2*L_z^2./(beta_z.^2.*(beta_z.^2+L_z^2+L_z));

            for iter_time = 1: numel(time)

                time_loc = time(iter_time);
                a_x = coeff_x.*exp(-(beta_x.^2*D*time_loc/half_thickness_x^2));    
                %a_y = coeff_y.*exp(-(beta_y.^2*D*time_loc/half_thickness_y^2));
                %a_z = coeff_z.*exp(-(beta_z.^2*D*time_loc/half_thickness_z^2));

                out_sigma_n(iter_time) = 1-sum(a_x).^3;
            end
        end
        
        % sigma_n for box with different side lengths
        function out_sigma_n = sigma_n_det_unsymm(obj, time)

            k = obj.k_ref;
            D = obj.D_ref;

            half_thickness_x = obj.half_thickness_x_ref;
            half_thickness_y = obj.half_thickness_y_ref;
            half_thickness_z = obj.half_thickness_z_ref;

            L_x = half_thickness_x*k/D;
            L_y = half_thickness_y*k/D;
            L_z = half_thickness_z*k/D;

            beta_x = obj.diffbeta(L_x);
            beta_y = obj.diffbeta(L_y);
            beta_z = obj.diffbeta(L_z);

            out_sigma_n = zeros(size(time));
            coeff_x = 2*L_x^2./(beta_x.^2.*(beta_x.^2+L_x^2+L_x));
            coeff_y = 2*L_y^2./(beta_y.^2.*(beta_y.^2+L_y^2+L_y));
            coeff_z = 2*L_z^2./(beta_z.^2.*(beta_z.^2+L_z^2+L_z));

            for iter_time = 1: numel(time)

                time_loc = time(iter_time);
                a_x = coeff_x.*exp(-(beta_x.^2*D*time_loc/half_thickness_x^2));    
                a_y = coeff_y.*exp(-(beta_y.^2*D*time_loc/half_thickness_y^2));
                a_z = coeff_z.*exp(-(beta_z.^2*D*time_loc/half_thickness_z^2));

                out_sigma_n(iter_time) = 1-sum(a_x)*sum(a_y)*sum(a_z);
            end
        end
        
        % sigma_n for box with 1D
        function out_sigma_n = sigma_n_det_1D(obj, time)

            k = obj.k_ref;
            D = obj.D_ref;

            half_thickness_x = obj.half_thickness_x_ref;
            half_thickness_y = obj.half_thickness_y_ref;
            half_thickness_z = obj.half_thickness_z_ref;
            
            half_thickness_x = min([half_thickness_x; half_thickness_y; half_thickness_z]);

            L_x = half_thickness_x*k/D;
            
            beta_x = obj.diffbeta(L_x);

            out_sigma_n = zeros(size(time));
            coeff_x = 2*L_x^2./(beta_x.^2.*(beta_x.^2+L_x^2+L_x));

            for iter_time = 1: numel(time)

                time_loc = time(iter_time);
                a_x = coeff_x.*exp(-(beta_x.^2*D*time_loc/half_thickness_x^2));

                out_sigma_n(iter_time) = 1-sum(a_x);
            end
        end
        
        % sigma_n for box with 2D        
        function out_sigma_n = sigma_n_det_2D(obj, time)

            k = obj.k_ref;
            D = obj.D_ref;

            half_thickness_x = obj.half_thickness_x_ref;
            half_thickness_y = obj.half_thickness_y_ref;
            half_thickness_z = obj.half_thickness_z_ref;
            
            dummy_half_thickness = [half_thickness_x; half_thickness_y; half_thickness_z];
            dummy_location = find(dummy_half_thickness<Inf);
            dummy_half_thickness = dummy_half_thickness(dummy_location);
            half_thickness_x = dummy_half_thickness(1);
            half_thickness_y = dummy_half_thickness(2);
            
            L_x = half_thickness_x*k/D;
            L_y = half_thickness_y*k/D;

            beta_x = obj.diffbeta(L_x);
            beta_y = obj.diffbeta(L_y);

            out_sigma_n = zeros(size(time));
            coeff_x = 2*L_x^2./(beta_x.^2.*(beta_x.^2+L_x^2+L_x));
            coeff_y = 2*L_y^2./(beta_y.^2.*(beta_y.^2+L_y^2+L_y));

            for iter_time = 1: numel(time)

                time_loc = time(iter_time);
                a_x = coeff_x.*exp(-(beta_x.^2*D*time_loc/half_thickness_x^2));    
                a_y = coeff_y.*exp(-(beta_y.^2*D*time_loc/half_thickness_y^2));

                out_sigma_n(iter_time) = 1-sum(a_x)*sum(a_y);
            end
        end
        
        % will add sphere
        
        function out_error_sigma = error_fct(obj, input_sigma_n_det)
            % we can modify the error structure later on
            out_error_sigma = obj.standard_dev*randn(size(input_sigma_n_det));
        end
        
    end
    
    % static functions
    methods (Static = true)
        function beta = diffbeta(l)


            % method of Fischer and Hertz. Solid State Ionics 2012

            M = 80;
            beta = zeros(M,1);
            if(~isnan(l))
                %options = optimset('Display', 'notify', 'FunValCheck', 'on', 'TolFun', 1e-10, 'TolX', 1E-12, 'MaxIter', 200, 'MaxFunEvals', 1000);
                options = optimset('FunValCheck', 'on', 'TolFun', 1e-17, 'TolX', 1E-20, 'MaxIter', 10000, 'MaxFunEvals', 100000);
                
                if l<1e-10
                    beta = [0:M-1]*pi;
                    beta(1) = fzero(@(x)(x*tan(x))-l,[0,(pi/2)-1e5*eps],options);
                elseif l>1e13
                    beta = (2*[0:M-1]+1)*(pi/2);
                else
                    
                    continue_solve_zeros = 1.;
                    
                    for n = 0:M-1
                        
                        if continue_solve_zeros
                            beta_local = fzero(@(x)(x*tan(x))-l,[(n*pi),(2*n+1)*(pi/2)-100*eps],options);
                        else
                            beta_local = n*pi;
                        end
                        
                        if abs(beta_local-n*pi)<1E-20
                            continue_solve_zeros = 0.0;
                        end
                        
                        beta(n+1) = beta_local;
                    end
                end
                
            end
            
        end;
        
    end
    
end   



