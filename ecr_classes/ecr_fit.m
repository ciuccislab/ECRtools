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


classdef ecr_fit < ecr_sens
    
    properties
        
        meas_time % measurement time
        meas_sigma_n % measurement sigma
        
    end
    
    methods
        
        function [theta_out] = log10_fit(obj)
            
            fun = @(log10_theta) dist_sigma_n_sigma_n_meas_log(obj, log10_theta);
            
            lb_log10_theta = [-7; -7];
            ub_log10_theta = [7; 7];
            
            %             % this works only with v2.0 of the OPTI Toolbox
            %             fprintf('Pre-fitting data using nlopt \n');
            %             Opt = opti('obj',fun,'bounds', lb_log10_theta, ub_log10_theta);
            
%             %nomad
%             fprintf('Pre-fitting data using nomad \n');
%             prob = optiprob('fun',fun,'bounds',lb_log10_theta,ub_log10_theta);
%             opts = optiset('solver','nomad');
%             Opt = opti(prob,opts);
            
            % pswarm
            fprintf('Pre-fitting data using pswarm \n');
            prob = optiprob('fun',fun,'bounds',lb_log10_theta,ub_log10_theta);
            opts = optiset('solver','pswarm');
            Opt = opti(prob,opts);
            
            log10_theta_0 = [0.0; 0.0];
            
            time_dummy = cputime;
            [log10_theta_out, dist_opt_nlopt, exitflag, info] = solve(Opt, log10_theta_0);
            time_dummy = cputime-time_dummy;
            fprintf('(time = %f s)\n', time_dummy);
            theta_out = 10.^log10_theta_out;
            
        end
        
        
        function [theta_out] = fit(obj, solver_used)
            
            fun = @(theta)dist_sigma_n_sigma_n_meas(obj, theta);
            
            lb_theta = [1E-3; 1E-3];
            ub_theta = [1E3; 1E3];
            
            fprintf('Fitting data using %s ', solver_used);
            theta_0 = [1.0; 1.0];

            switch solver_used
                case 'pswarm'
                    
                    prob = optiprob('fun',fun,'bounds',lb_theta,ub_theta);
                    opts = optiset('solver','pswarm','display','iter', 'maxiter', 20000);
                    Opt = opti(prob,opts);
                    
                case 'nomad'
                    
                    prob = optiprob('fun',fun,'bounds',lb_theta,ub_theta);
                    opts = optiset('solver','nomad');
                    Opt = opti(prob,opts);
                    
                case 'GN_DIRECT'
                    
                    prob = optiprob('fun',fun,'bounds',lb_theta,ub_theta);
                    opts = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'));
                    Opt = opti(prob,opts);
                    
                case 'nlopt'
                    
                    nopt = nloptset('algorithm','LN_BOBYQA', 'subtolafun', 1E-8, 'subtolrfun', 1E-6);            % nloptset functions much like optiset
                    opts = optiset('solver','nlopt','solverOpts',nopt, 'tolrfun', 1E-9, 'tolafun', 1E-9, 'tolint', 1E-8);  % specify solver + our NLOPT options above
                    Opt = opti('obj',fun,'bounds', lb_theta, ub_theta,'ndec',2,'options',opts);
                    theta_0 = [1.0; 1.0] + 0.01*randn(2,1);
                    
            end
                     
            time_dummy = cputime;
            [theta_out, dist_opt_nlopt] = solve(Opt, theta_0);
            time_dummy = cputime-time_dummy;
            fprintf(' (time = %f s)\n', time_dummy);
            
        end
   
        
        function theta_vec = synthetic_exp(obj, solver_used, N_exp)
            
            time_exp = obj.meas_time;
            base_sigma_n_det = sigma_n_det(obj, time_exp);
            
            theta_vec = zeros(2, N_exp);
            
            for iter_synth_exp = 1: N_exp
                
                loc_error_sigma = error_fct(obj, base_sigma_n_det);
                loc_meas_sigma_n = base_sigma_n_det+ loc_error_sigma;
                obj.meas_sigma_n = loc_meas_sigma_n;
                
                fun = @(theta)dist_sigma_n_sigma_n_meas(obj, theta);
                lb_theta = [1E-1; 1E-1];
                ub_theta = [1E1; 1E1];
                theta_0 = [1;1];
                
                fprintf('Synthetic measurement %u/%u (fit done with %s) ', iter_synth_exp, N_exp, solver_used);
                
                
                switch solver_used
                    case 'pswarm'
                        
                        prob = optiprob('fun',fun,'bounds',lb_theta,ub_theta);
                        opts = optiset('solver','pswarm');
                        Opt = opti(prob,opts);
                        
                    case 'nomad'
                        
                        prob = optiprob('fun',fun,'bounds',lb_theta,ub_theta);
                        opts = optiset('solver','nomad');
                        Opt = opti(prob,opts);
                        
                    case 'GN_DIRECT'
                        
                        prob = optiprob('fun',fun,'bounds',lb_theta,ub_theta);
                        opts = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'));
                        Opt = opti(prob,opts);
                        
                    case 'nlopt'
                        
                        nopt = nloptset('algorithm','LN_BOBYQA', 'subtolafun', 1E-8, 'subtolrfun', 1E-7);            % nloptset functions much like optiset
                        opts = optiset('solver','nlopt','solverOpts',nopt, 'tolrfun', 1E-9, 'tolafun', 1E-9, 'tolint', 1E-8);  % specify solver + our NLOPT options above
                        Opt = opti('obj',fun,'bounds', lb_theta, ub_theta,'ndec',2,'options',opts);
                        
                end
                
                
                time_dummy = cputime;
                [loc_theta,fval,exitflag,info] = solve(Opt,theta_0);
                time_dummy = cputime-time_dummy;
                fprintf('(time = %f s)\n', time_dummy);
                
                theta_vec(:, iter_synth_exp) = loc_theta;
                
                
            end
            
        end
        
        
        function output_sigma_n_det = output_sigma_n_det(obj, theta)
            
            loc_meas_time = obj.meas_time;
            output_sigma_n_det = sigma_n_det_theta(obj, theta, loc_meas_time);
            
        end
        
        function output_std = compute_std_meas(obj)
            
            theta = [1; 1];
            sigma_n_det = output_sigma_n_det(obj, theta);
            loc_meas_sigma_n = obj.meas_sigma_n;
            
            computed_variance = 1/numel(sigma_n_det)*norm(sigma_n_det-loc_meas_sigma_n)^2;
            computed_std = sqrt(computed_variance);
            
            obj.standard_dev = computed_std;
            output_std = computed_std;
            
            
        end
        
        
        
    end
    
    methods (Access = protected)
        
        function out_dist = dist_sigma_n_sigma_n_meas_log(obj, log10_theta)
            
            loc_meas_time = obj.meas_time;
            loc_meas_sigma_n = obj.meas_sigma_n;
            
            theta = 10.^log10_theta;
            sigma_n_est = sigma_n_det_theta(obj, theta, loc_meas_time);
            
            out_dist = norm(sigma_n_est-loc_meas_sigma_n);
            out_dist = out_dist^2;
            
        end
        
        
        function out_dist = dist_sigma_n_sigma_n_meas(obj, theta)
            
            loc_meas_time = obj.meas_time;
            loc_meas_sigma_n = obj.meas_sigma_n;
            
            sigma_n_est = sigma_n_det_theta(obj, theta, loc_meas_time);
            
            out_dist = norm(sigma_n_est-loc_meas_sigma_n);
            out_dist = out_dist^2;
            
        end
        
    end
end



