classdef spatial_express_superphase < bt.model.template_spatial
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
	end

	methods
		function self = spatial_express_superphase(n_modes) % Constructor
			self.name = 'spatial_express_superphase';
			self.param_names = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','Alpha','Beta','t_0 ','t_0 amp','t_0 phase','G_{ee} amp','G_{ee} phase'};
			self.n_params = length(self.param_names);
			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.1.*[
				10    	% Gee
				10    	% Gei
				10    	% Gese
				20    	% Gesre
				20    	% Gsrs
				40    	% Alpha
				100   	% Beta
				0.02 	% t0
				0.01 	% t0_amp
				0.1  	% t0_phase
				2    	% gee_amp
				0.1  	% gee_phase
			].';
			self.limits = [
				eps       	20     	  %Gee
				-20     	-eps      %Gei
				eps       	40     	  %Gese
				-40     	-eps      %Gesre
				-14     	-eps      %Gsrs
				10      	100    	  %Alpha
				100     	800    	  %Beta
				0.050  		0.1400 	  %t0
				0           0.05      %t0_amp
				0 			0.5		  % t0_phase
				0 			5		  % gee_amp
				0 			0.5		  % gee_phase
			].';
			assert(numel(self.limits) == self.n_params*2)
			self.p = model.params();
			self.p.disable_set = true;
		end

		function w = get_weights(self,target_f) % Default weighting function
			w = get_weights@bt.model.template(self,target_f);
			w(target_f > 25) = 0;
		end
		
		function valid = validate_params(self,pars) % Check if parameters are valid
			% This implements a basic range check. Derived classes should impose their own checks as well
			% Probably simplest to overload this function
			valid = ~(pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || pars(7)/pars(6) > 20 || any(pars > self.limits(2,:) | pars < self.limits(1,:) ) );

			% Also validate the cosine parameters
			valid = valid && self.check_cosine_validity(pars,8,9) && self.check_cosine_validity(pars,1,11);
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
			self.p.emg_a = 0;
			self.p.apply_variation('t0','cosine_phased',pars(9),pars(10));
			self.p.apply_variation('g5_ee','cosine_phased',pars(11),pars(12));
		end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,mean(target_P,2),0); % Fit the first (or only) spectrum given
			a =  db_data.gab(idx,:);
			p = model.params(db_data.iswake(idx));
			initial_values =  [a(1) a(2) a(3)*a(4) a(3)*a(5)*a(7) a(5)*a(8) p.alpha(1) p.beta(1) p.t0 0 0 0 0];
			prior_pp = self.uniform_priors();
		end

		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
		end
		
	end
end



