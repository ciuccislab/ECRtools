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


classdef ecr_sens < ecr_base
    
    % methods/functions in the class accessible to the user
    methods
        
            function out_grad_sigma_n = grad_sigma_n(obj, time)
               
                delta = 1E-3;
                
                theta_k = 1.0 + delta; theta_D = 1.0;
                theta = [theta_k; theta_D];
                sigma_n_delta_theta_k = sigma_n_det_theta(obj, theta, time);
                
                theta_k = 1.0 - delta; theta_D = 1.0;
                theta = [theta_k; theta_D];
                sigma_n_m_delta_theta_k = sigma_n_det_theta(obj, theta, time);
                
                theta_k = 1.0; theta_D = 1.0 + delta;
                theta = [theta_k; theta_D];
                sigma_n_delta_theta_D = sigma_n_det_theta(obj, theta, time);
                
                theta_k = 1.0; theta_D = 1.0 - delta;
                theta = [theta_k; theta_D];
                sigma_n_m_delta_theta_D = sigma_n_det_theta(obj, theta, time);
                
                grad_sigma_n_theta_k = (sigma_n_delta_theta_k - sigma_n_m_delta_theta_k)/(2.*delta);
                grad_sigma_n_theta_D = (sigma_n_delta_theta_D - sigma_n_m_delta_theta_D)/(2.*delta);

                out_grad_sigma_n = [grad_sigma_n_theta_k, grad_sigma_n_theta_D];
                
            end
            
            function matrix_V = cov_sigma_n(obj, time)
                
                
                out_grad_sigma_n = grad_sigma_n(obj, time);
                
                grad_sigma_n_theta_k = out_grad_sigma_n(:, 1);
                grad_sigma_n_theta_D = out_grad_sigma_n(:, 2);

                A_theta = zeros(2);
                
                A_theta(1,1) = grad_sigma_n_theta_k'*grad_sigma_n_theta_k;
                A_theta(1,2) = grad_sigma_n_theta_k'*grad_sigma_n_theta_D;
                A_theta(2,1) = grad_sigma_n_theta_D'*grad_sigma_n_theta_k;
                A_theta(2,2) = grad_sigma_n_theta_D'*grad_sigma_n_theta_D;

                if condest(A_theta)>1E25
                    matrix_V = 10^20*eye(2);
                else
                    matrix_V = (obj.standard_dev)^2*inv(A_theta);
                end                
  
            end
            
            function out_det = det_cov_sigma_n(obj, time)
                
                matrix_V = cov_sigma_n(obj, time);
                out_det = det(matrix_V);
                
            end
            
            
            function out_trace = trace_cov_sigma_n(obj, time)
                
                matrix_V = cov_sigma_n(obj, time);
                out_trace = trace(matrix_V);
                
            end
            
            function out_cond = cond_cov_sigma_n(obj, time)
                
                matrix_V = cov_sigma_n(obj, time);
                out_cond = condest(matrix_V);
                
            end
            
            function out_max_eig = max_eig_cov_sigma_n(obj, time)
                
                matrix_V = cov_sigma_n(obj, time);
                out_max_eig = eigs(matrix_V,1);
                
            end

            function data_boundary_theta  = boundary_cov_sigma_n(obj, time)

                matrix_V = cov_sigma_n(obj, time);
                inv_matrix_V = inv(matrix_V);

                [X_mat, lambdas] = eig(inv_matrix_V);
                theta = linspace(0, 2*pi, 101);
                x_prime = 3/sqrt(lambdas(1,1))*cos(theta);
                y_prime = 3/sqrt(lambdas(2,2))*sin(theta);

                xy_prime = [x_prime; y_prime];
                xy = X_mat*xy_prime;

                xy(1, :) = xy(1, :) + 1.0;
                xy(2, :) = xy(2, :) + 1.0;

                data_boundary_theta = xy;

            end
            
            function out_sigma_n = sigma_n_det_theta(obj, theta, time)
                
                half_thickness_x = obj.half_thickness_x_ref;
                half_thickness_y = obj.half_thickness_y_ref;
                half_thickness_z = obj.half_thickness_z_ref;
                
                theta_k = theta(1);
                theta_D = theta(2);
                
                if half_thickness_x == half_thickness_y && half_thickness_x == half_thickness_z
                    out_sigma_n = sigma_n_det_symm_theta(obj, theta_k, theta_D, time);
                else
                    out_sigma_n = sigma_n_det_unsymm_theta(obj, theta_k, theta_D, time);
                end
                
                
                if half_thickness_x == Inf || half_thickness_y == Inf || half_thickness_z == Inf
                    if half_thickness_x == Inf && half_thickness_y<Inf && half_thickness_z<Inf
                        out_sigma_n = sigma_n_det_2D_theta(obj, theta_k, theta_D, time);
                    elseif half_thickness_y == Inf && half_thickness_x<Inf && half_thickness_z<Inf
                        out_sigma_n = sigma_n_det_2D_theta(obj, theta_k, theta_D, time);
                    elseif half_thickness_z == Inf && half_thickness_x<Inf && half_thickness_y<Inf
                        out_sigma_n = sigma_n_det_2D_theta(obj, theta_k, theta_D, time);
                    else
                        out_sigma_n = sigma_n_det_1D_theta(obj, theta_k, theta_D, time);
                    end
                elseif half_thickness_x == half_thickness_y && half_thickness_x == half_thickness_z
                    out_sigma_n = sigma_n_det_symm_theta(obj, theta_k, theta_D, time);
                else
                    out_sigma_n = sigma_n_det_unsymm_theta(obj, theta_k, theta_D, time);
                end
            end

    end

    methods (Access = protected)
        % sigma_n for box with different side lengths
        function out_sigma_n = sigma_n_det_unsymm_theta(obj, theta_k, theta_D, time)

            k = theta_k*obj.k_ref;
            D = theta_D*obj.D_ref;

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
        
        function out_sigma_n = sigma_n_det_symm_theta(obj, theta_k, theta_D, time)

            k = theta_k*obj.k_ref;
            D = theta_D*obj.D_ref;

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
        
        % sigma_n for box with 1D
        function out_sigma_n = sigma_n_det_1D_theta(obj, theta_k, theta_D, time)

            k = obj.k_ref*theta_k;
            D = obj.D_ref*theta_D;

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
        function out_sigma_n = sigma_n_det_2D_theta(obj, theta_k, theta_D, time)

            k = obj.k_ref*theta_k;
            D = obj.D_ref*theta_D;

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

    end
    
    
end         
        

        