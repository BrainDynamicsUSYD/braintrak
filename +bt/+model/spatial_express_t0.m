classdef spatial_express_t0 < bt.model.template_spatial
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
	end

	methods
		function self = spatial_express_t0(n_modes) % Constructor
			self.name = 'spatial_express_t0';
			self.n_params = 10;
			self.param_names = {'Gee','Gei','Gese','Gesre','Gsrs','Alpha','Beta','t0','EMGa','t0_amp'};
			self.param_symbols = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta','t_0','A_{EMG}','t_0amp'};
			self.param_units = {'','','','','','s^{-1}','s^{-1}','ms','','ms'};

			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.05.*[2 2 2 2 2 5 10 0.02 0.01 0.05];
			self.limits = [
				eps       	20     	  %Gee
				-20     	-eps      %Gei
				eps       	40     	  %Gese
				-40     	-eps      %Gesre
				-14     	-eps      %Gsrs
				10      	100    	  %Alpha
				100     	800    	  %Beta
				0.0750  	0.1400 	  %t0
				0       	5      	  %EMGa
				-1       	1      	  %t0_amp
				].';

			self.p = model.params();
			self.p.disable_set = true;
		end


		function valid = validate_params(self,pars) % Check if parameters are valid
			% This implements a basic range check. Derived classes should impose their own checks as well
			% Probably simplest to overload this function
			valid = ~(pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || pars(7)/pars(6) > 20 || any(pars > self.limits(2,:) | pars < self.limits(1,:) ) );
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
			self.p.emg_a = pars(9);
			self.p.spatial_t0 = self.p.variation_function('cosine',self.p.t0,pars(10));
		end

		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,mean(target_P,2)); % Fit the first (or only) spectrum given
			a =  db_data.gab(idx,:);
			p = model.params(db_data.iswake(idx));
			initial_values =  [a(1) a(2) a(3)*a(4) a(3)*a(5)*a(7) a(5)*a(8) p.alpha(1) p.beta(1) p.t0 0 0];
			prior_pp = self.uniform_priors();
		end
		

		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
		end
		
	end
end



