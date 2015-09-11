classdef spatial_emg < mcmc.model.template_spatial
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
	end

	methods
		function self = spatial_emg(n_modes) % Constructor
			self.name = 'spatial_emg';
			self.param_names = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','Alpha','Beta','t_0 ','t_0 amp','t_0 phase','EMG_a','EMG_x amp','EMG_x phase','EMG_y amp','EMG_y phase'};
			self.param_symbols = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta','t_0','t_0_amp','t_0 phase','EMG_a','EMG_x amp','EMG_x phase','EMG_y amp','EMG_y phase'};
			self.param_units = {'','','','','','s^{-1}','s^{-1}','ms','ms','m','','','','',''};

			self.n_params = length(self.param_names);
			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.1.*[
				0.4    	% Gee
				0.4    	% Gei
				1    	% Gese
				1    	% Gesre
				0.2    	% Gsrs
				5    	% Alpha
				40   	% Beta
				0.005 	% t0
				0.01 	% t0_amp
				0.25  	% t0_phase
				0.0001 	% emg_a
				0.0001 	% emg_x_var
				0.25 	% emg_x_phase
				0.0001 	% emg_y_var
				0.25 	% emg_y_phase
			].';

			Lx = 0.5;
			head_d = mcmc.electrode_positions;
			offset = (Lx-head_d)/2;

			self.limits = [
				eps       	20     	  %Gee
				-20     	-eps      %Gei
				eps       	40     	  %Gese
				-40     	-eps      %Gesre
				-14     	-eps      %Gsrs
				10      	100    	  %Alpha
				100     	800    	  %Beta
				0.050  		0.1400 	  %t0
				-0.05           0.05      %t0_amp
				0+offset 	0.5-offset		  % t0_phase
				0 			5		  % emg_a
				-1 			1		  % emg_x_var
				0+offset 	0.5-offset  % emg_x_phase
				-1 			1		  % emg_y_var
				0+offset 	0.5-offset	% emg_y_phase
			].';

			self.emg_scaling = 1e-9;

			assert(numel(self.limits) == self.n_params*2)
			self.p = model.params();
			self.p.phin = 1e-5;
			self.p.disable_set = true;
		end

		% function w = get_weights(self,target_f) % Default weighting function
		% 	w = get_weights@mcmc.model.template(self,target_f);
		% 	w(target_f > 25) = 0;
		% end
		
		function valid = validate_params(self,pars) % Check if parameters are valid
			% This implements a basic range check. Derived classes should impose their own checks as well
			% Probably simplest to overload this function
			valid = ~(pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || pars(7)/pars(6) > 20 || any(pars > self.limits(2,:) | pars < self.limits(1,:) ) );
			
			% Also validate the cosine parameters
			valid = valid && self.check_cosine_validity(pars,8,9) && self.check_cosine_validity(pars,11,12) && self.check_cosine_validity(pars,11,14);
			%valid
		end

		function params = params_from_p(self,p) % Map point to parameters
			params = [p.gabcd p.alpha(1) p.beta(1) p.t0 NaN NaN p.emg_a/self.emg_scaling NaN NaN NaN NaN];
		end

		function set_params(self,pars)
			% This function is what decides what the spatial variations are actually going to be
			%e.g. p.apply_variation('cosine','t0',pars(10));
			self.p.gabcd = pars(1:5);
			self.p.alpha(:) = pars(6);
			self.p.beta(:) = pars(7);
			self.p.t0 = pars(8);
			self.p.taues = self.p.t0/2;
			self.p.tause = self.p.t0/2;
			self.p.emg_a = self.emg_scaling*pars(11);
			self.p.apply_variation('t0','cosine_phased',pars(9),pars(10));
			self.p.apply_variation('emg','cosine_2d_phased',self.emg_scaling*pars(12),pars(13),self.emg_scaling*pars(14),pars(15));
		end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,mean(target_P,2),0); % Fit the first (or only) spectrum given
			a =  db_data.gab(idx,:);
			p = model.params(db_data.iswake(idx));
			initial_values =  [a(1) a(2) a(3)*a(4) a(3)*a(5)*a(7) a(5)*a(8) p.alpha(1) p.beta(1) p.t0 0 0.25 0 0 0.25 0 0.25];
			prior_pp = self.uniform_priors();
		end

		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
		end
		
	end
end



