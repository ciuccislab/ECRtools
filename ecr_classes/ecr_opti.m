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


classdef ecr_opti < ecr_fit
    
    properties
    
        
    end
    
    methods
            
        function [time_exp_opt, thickness_opt, goal_opt] = optimize_exp_Nps(obj, time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, dimension)
            
            xi_0 = [1; 1];
            time_exp_ref = (time_exp_bound(1)+time_exp_bound(2))/2.0;
            half_thickness_ref = (half_thickness_bound(1)+half_thickness_bound(2))/2.0;
            lb_xi = [time_exp_bound(1)/time_exp_ref; half_thickness_bound(1)/half_thickness_ref];
            ub_xi = [time_exp_bound(2)/time_exp_ref; half_thickness_bound(2)/half_thickness_ref];
            
            xi_diml_ref = [time_exp_ref; half_thickness_ref];
            
            switch criterion
                case 'trace'
                    fun = @(xi) opt_trace_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);
                
                case 'max_eig' 
                    fun = @(xi) opt_max_eig_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);
                
                case 'det'
                    fun = @(xi) opt_det_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);
                
                case 'cond'
                    fun = @(xi) opt_cond_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);
            end
            
            prob = optiprob('fun', fun, 'bounds', lb_xi, ub_xi);

            % nlopt
           
            opts = optiset('solver', 'pswarm'); 
            Opt = opti(prob,opts);
            [xi_opt, goal_opt, exitflag, info] = solve(Opt, xi_0);
            
            time_exp_opt = xi_opt(1)*time_exp_ref;
            thickness_opt = xi_opt(2)*half_thickness_ref;
            goal_opt = 10^goal_opt;
            
        end
        
        function [time_exp_opt, thickness_opt, goal_opt] = optimize_exp_Ntot(obj, time_exp_bound, half_thickness_bound, N_tot, criterion, dimension)
            
            xi_0 = [1; 1];
            time_exp_ref = (time_exp_bound(1)+time_exp_bound(2))/2.0;
            half_thickness_ref = (half_thickness_bound(1)+half_thickness_bound(2))/2.0;
            lb_xi = [time_exp_bound(1)/time_exp_ref; half_thickness_bound(1)/half_thickness_ref];
            ub_xi = [time_exp_bound(2)/time_exp_ref; half_thickness_bound(2)/half_thickness_ref];
            
            xi_diml_ref = [time_exp_ref; half_thickness_ref];
            
            switch criterion
                case 'trace'
                    fun = @(xi) opt_trace_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);
                
                case 'max_eig' 
                    fun = @(xi) opt_max_eig_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);
                
                case 'det'
                    fun = @(xi) opt_det_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);        
                
                case 'cond'
                    fun = @(xi) opt_cond_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);
            end
            
            prob = optiprob('fun', fun, 'bounds', lb_xi, ub_xi);

            % nlopt
           
            opts = optiset('solver', 'pswarm'); 
            Opt = opti(prob,opts);
            [xi_opt, goal_opt, exitflag, info] = solve(Opt, xi_0);
            
            time_exp_opt = xi_opt(1)*time_exp_ref;
            thickness_opt = xi_opt(2)*half_thickness_ref;
            goal_opt = 10^goal_opt;
            
        end
        
    %%%%%%%%%%% Robust    
        function [time_exp_opt, thickness_opt, goal_opt] = optimize_avg_robust_Ntot(obj, time_exp_bound, half_thickness_bound, N_tot, criterion, dimension, k_int, D_int)
            
            xi_0 = [1; 1];
            time_exp_ref = (time_exp_bound(1)+time_exp_bound(2))/2.0;
            half_thickness_ref = (half_thickness_bound(1)+half_thickness_bound(2))/2.0;
            lb_xi = [time_exp_bound(1)/time_exp_ref; half_thickness_bound(1)/half_thickness_ref];
            ub_xi = [time_exp_bound(2)/time_exp_ref; half_thickness_bound(2)/half_thickness_ref];
            
            xi_diml_ref = [time_exp_ref; half_thickness_ref];
            
            fun = @(xi) avg_robust_integral_phi_Ntot(obj, xi, xi_diml_ref, N_tot, criterion, dimension, k_int, D_int);
            
            prob = optiprob('fun', fun, 'bounds', lb_xi, ub_xi);

            % nlopt
           
            opts = optiset('solver', 'pswarm'); 
            Opt = opti(prob,opts);
            tic
            [xi_opt, goal_opt, exitflag, info] = solve(Opt, xi_0);
            toc
            time_exp_opt = xi_opt(1)*time_exp_ref;
            thickness_opt = xi_opt(2)*half_thickness_ref;
            
            goal_opt = 10^goal_opt;
            
        end
         
        function [time_exp_opt, thickness_opt, goal_opt] = optimize_avg_robust_Nps(obj, time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, dimension, k_int, D_int)
            
            xi_0 = [1; 1];
            time_exp_ref = (time_exp_bound(1)+time_exp_bound(2))/2.0;
            half_thickness_ref = (half_thickness_bound(1)+half_thickness_bound(2))/2.0;
            lb_xi = [time_exp_bound(1)/time_exp_ref; half_thickness_bound(1)/half_thickness_ref];
            ub_xi = [time_exp_bound(2)/time_exp_ref; half_thickness_bound(2)/half_thickness_ref];
            
            xi_diml_ref = [time_exp_ref; half_thickness_ref];
            

            fun = @(xi) avg_robust_integral_phi_Nps(obj, xi, xi_diml_ref, N_meas_per_second, criterion, dimension, k_int, D_int) 
            
            prob = optiprob('fun', fun, 'bounds', lb_xi, ub_xi);

            % nlopt
           
            opts = optiset('solver', 'pswarm'); 
            Opt = opti(prob,opts);
            tic
            [xi_opt, goal_opt, exitflag, info] = solve(Opt, xi_0);
            toc
            time_exp_opt = xi_opt(1)*time_exp_ref;
            thickness_opt = xi_opt(2)*half_thickness_ref;
            
            goal_opt = 10^goal_opt;
            
        end
        
        function [time_exp_opt, thickness_opt, goal_opt] = optimize_minmax_robust_Ntot(obj, time_exp_bound, half_thickness_bound, N_tot, criterion, dimension, k_int, D_int)
            
            xi_0 = [1; 1];
            time_exp_ref = (time_exp_bound(1)+time_exp_bound(2))/2.0;
            half_thickness_ref = (half_thickness_bound(1)+half_thickness_bound(2))/2.0;
            lb_xi = [time_exp_bound(1)/time_exp_ref; half_thickness_bound(1)/half_thickness_ref];
            ub_xi = [time_exp_bound(2)/time_exp_ref; half_thickness_bound(2)/half_thickness_ref];
            
            xi_diml_ref = [time_exp_ref; half_thickness_ref];
            

            fun = @(xi) max_robust_phi_Ntot(obj, xi, xi_diml_ref, N_tot, criterion, dimension, k_int, D_int);
            
            prob = optiprob('fun', fun, 'bounds', lb_xi, ub_xi);

            % nlopt
           
            opts = optiset('solver', 'pswarm'); 
            %opts = optiset('solver', 'nlopt'); 
            
            Opt = opti(prob,opts);
            tic
            [xi_opt, goal_opt, exitflag, info] = solve(Opt, xi_0);
            toc
            time_exp_opt = xi_opt(1)*time_exp_ref;
            thickness_opt = xi_opt(2)*half_thickness_ref;
            
            goal_opt = 10^goal_opt;
            
        end
       
        function [time_exp_opt, thickness_opt, goal_opt] = optimize_minmax_robust_Nps(obj, time_exp_bound, half_thickness_bound, N_meas_per_second, criterion, dimension, k_int, D_int)
            
            xi_0 = [1; 1];
            time_exp_ref = (time_exp_bound(1)+time_exp_bound(2))/2.0;
            half_thickness_ref = (half_thickness_bound(1)+half_thickness_bound(2))/2.0;
            lb_xi = [time_exp_bound(1)/time_exp_ref; half_thickness_bound(1)/half_thickness_ref];
            ub_xi = [time_exp_bound(2)/time_exp_ref; half_thickness_bound(2)/half_thickness_ref];
            
            xi_diml_ref = [time_exp_ref; half_thickness_ref];
            

            fun = @(xi) max_robust_phi_Nps(obj, xi, xi_diml_ref, N_meas_per_second, criterion, dimension, k_int, D_int);
            
            prob = optiprob('fun', fun, 'bounds', lb_xi, ub_xi);

            % nlopt
           
            %opts = optiset('solver', 'pswarm'); 
            opts = optiset('solver', 'nlopt'); 
            
            Opt = opti(prob,opts);
            tic
            [xi_opt, goal_opt, exitflag, info] = solve(Opt, xi_0);
            toc
            time_exp_opt = xi_opt(1)*time_exp_ref;
            thickness_opt = xi_opt(2)*half_thickness_ref;
            
            goal_opt = 10^goal_opt;
            
        end
            
    end
    
    methods (Access = protected)
        
        function out_det = opt_det_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);

                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_det = log10(det_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_max_eig = opt_max_eig_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_max_eig = log10(max_eig_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_trace = opt_trace_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_trace = log10(trace_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_cond = opt_cond_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_cond = log10(cond_cov_sigma_n(obj, time_exp_local));
        end
            
        function out_det = opt_det_Ntot(obj, xi, xi_diml_ref, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_det = log10(det_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_max_eig = opt_max_eig_Ntot(obj, xi, xi_diml_ref, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_max_eig = log10(max_eig_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_trace = opt_trace_Ntot(obj, xi, xi_diml_ref, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_trace = log10(trace_cov_sigma_n(obj, time_exp_local));
        end
         
        function out_cond = opt_cond_Ntot(obj, xi, xi_diml_ref, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                out_cond = log10(cond_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_det = max_det_Nps(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_det = -log10(det_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_max_eig = max_max_eig_Nps(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_max_eig = -log10(max_eig_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_trace = max_trace_Nps(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_trace = -log10(trace_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_cond = max_cond_Nps(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = [0:1/N_meas_per_second:time_exp_max]';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_cond = -log10(cond_cov_sigma_n(obj, time_exp_local));
        end
            
        function out_det = max_det_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_det = -log10(det_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_max_eig = max_max_eig_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_max_eig = -log10(max_eig_cov_sigma_n(obj, time_exp_local));
        end
        
        function out_trace = max_trace_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_trace = -log10(trace_cov_sigma_n(obj, time_exp_local));
        end
         
        function out_cond = max_cond_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension)    
                
                time_exp_max = floor(xi_diml_ref(1)*xi(1));
                time_exp_local = linspace(0,time_exp_max, N_tot)';
                
                half_thickness_loc = xi_diml_ref(2)*xi(2);
                
                switch dimension
                    case '1D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = Inf;
                        obj.half_thickness_z_ref = Inf;
                    case '2D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = Inf;
                    case '3D'
                        obj.half_thickness_x_ref = half_thickness_loc;
                        obj.half_thickness_y_ref = half_thickness_loc;
                        obj.half_thickness_z_ref = half_thickness_loc;
                end
                
                obj.k_ref = obj.k_ref*theta(1);
                obj.D_ref = obj.D_ref*theta(2);
                
                out_cond = -log10(cond_cov_sigma_n(obj, time_exp_local));
        end
        
        function integral_out = avg_robust_integral_phi_Ntot(obj, xi, xi_diml_ref, N_tot, criterion, dimension, k_int, D_int) 

            [k_int_mat, D_int_mat] = meshgrid(k_int, D_int);
            coeff = (max(k_int)-min(k_int))*(max(D_int)-min(D_int));

            phi_mat = zeros(size(k_int_mat));

            for iter_k = 1: numel(k_int)
    
                for iter_D = 1: numel(D_int)
        
                    k_ref = k_int_mat(iter_D, iter_k);
                    D_ref = D_int_mat(iter_D, iter_k);
                    
                    obj.k_ref = k_ref;
                    obj.D_ref = D_ref;
        
        
                    switch criterion
                        case 'trace'
                            phi_loc = opt_trace_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);
                
                        case 'max_eig' 
                            phi_loc = opt_max_eig_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);
                
                        case 'det'
                            phi_loc = opt_det_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);        
                
                        case 'cond'
                            phi_loc = opt_cond_Ntot(obj, xi, xi_diml_ref, N_tot, dimension);
                    end
                                                                       
                    phi_mat(iter_D, iter_k) = 10^phi_loc;
            
                end
            end

            integral_out = trapz(k_int, trapz(D_int, phi_mat))/coeff;
            integral_out = log10(integral_out);
            %fprintf('Integral =%f\n', integral_out);

                    
        end
                        
        function integral_out = avg_robust_integral_phi_Nps(obj, xi, xi_diml_ref, N_meas_per_second, criterion, dimension, k_int, D_int) 

            [k_int_mat, D_int_mat] = meshgrid(k_int, D_int);
            coeff = (max(k_int)-min(k_int))*(max(D_int)-min(D_int));

            phi_mat = zeros(size(k_int_mat));

            for iter_k = 1: numel(k_int)
    
                for iter_D = 1: numel(D_int)
        
                    k_ref = k_int_mat(iter_D, iter_k);
                    D_ref = D_int_mat(iter_D, iter_k);
                    
                    obj.k_ref = k_ref;
                    obj.D_ref = D_ref;
        
        
                    switch criterion
                        case 'trace'
                            phi_loc = opt_trace_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);
                
                        case 'max_eig' 
                            phi_loc = opt_max_eig_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);
                
                        case 'det'
                            phi_loc = opt_det_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);        
                
                        case 'cond'
                            phi_loc = opt_cond_Nps(obj, xi, xi_diml_ref, N_meas_per_second, dimension);
                    end
                                                                       
                    phi_mat(iter_D, iter_k) = 10^phi_loc;
            
                end
            end

            integral_out = trapz(k_int, trapz(D_int, phi_mat))/coeff;
            integral_out = log10(integral_out);
            fprintf('Integral =%f\n', integral_out);
        end
        
        function out_max = max_robust_phi_Ntot(obj, xi, xi_diml_ref, N_tot, criterion, dimension, k_int, D_int) 
                    
            here_k_ref = 1/2*(max(k_int)+min(k_int));
            here_D_ref = 1/2*(max(D_int)+min(D_int));
            
            obj.k_ref = here_k_ref;
            obj.D_ref = here_D_ref;
            
            switch criterion
                case 'trace'
                    fun = @(theta) max_trace_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension);
                    
                case 'max_eig'
                    fun = @(theta) max_max_eig_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension);
                    
                case 'det'
                    fun = @(theta) max_det_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension);
                    
                case 'cond'
                    fun = @(theta) max_cond_Ntot(obj, xi, xi_diml_ref, theta, N_tot, dimension);
            end
            
            lb_theta = [min(k_int)/here_k_ref; min(D_int)/here_D_ref];
            ub_theta = [max(k_int)/here_k_ref; max(D_int)/here_D_ref];
            theta_0 = [1;1];
            
            opts = optiset('solver', 'nomad');
            prob = optiprob('fun', fun, 'bounds', lb_theta, ub_theta);
            Opt = opti(prob,opts);
            
            [theta_opt, goal_opt, exitflag, info] = solve(Opt, theta_0);
            
            out_max = -goal_opt;
                    
        end

       function out_max = max_robust_phi_Nps(obj, xi, xi_diml_ref, N_meas_per_second, criterion, dimension, k_int, D_int) 
                    
            here_k_ref = 1/2*(max(k_int)+min(k_int));
            here_D_ref = 1/2*(max(D_int)+min(D_int));
            
            obj.k_ref = here_k_ref;
            obj.D_ref = here_D_ref;

            switch criterion
                case 'trace'
                    fun = @(theta) max_trace_Nps(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension);
                    
                case 'max_eig'
                    fun = @(theta) max_max_eig_Nps(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension);
                    
                case 'det'
                    fun = @(theta) max_det_Nps(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension);
                    
                case 'cond'
                    fun = @(theta) max_cond_Ntot(obj, xi, xi_diml_ref, theta, N_meas_per_second, dimension);
            end
            
            lb_theta = [min(k_int)/here_k_ref; min(D_int)/here_D_ref];
            ub_theta = [max(k_int)/here_k_ref; max(D_int)/here_D_ref];
            theta_0 = [1;1];
            
            opts = optiset('solver', 'nomad'); 
            prob = optiprob('fun', fun, 'bounds', lb_theta, ub_theta);
            Opt = opti(prob,opts);
            
            [theta_opt, goal_opt, exitflag, info] = solve(Opt, theta_0);
            
            out_max = -goal_opt;
                    
        end
        
    end
    
    
end 
        
        

        