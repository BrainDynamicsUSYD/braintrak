classdef full_b_ratio < bt.model.full
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit

	properties
	end


	methods
		function self = full_b_ratio() % Constructor
			self.name = 'model_full_b_ratio';
			self.n_params = 9;
			self.param_names = {'Gee','Gei','Gese','Gesre','Gsrs','Alpha','Beta/Alpha','t0','EMGa'};
			self.param_symbols = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta/\alpha','t_0','A_{EMG}'};
			self.param_units = {'','','','','','s^{-1}','','ms',''};

			self.n_fitted = 9;
			self.skip_fit = zeros(1,self.n_fitted);

			self.initial_step_size = [0.4  0.4  1  1  0.2  5  0.2  0.005  0.0001];
		   	self.limits = [ eps  -20        eps     -40      -14      10      1    0.075   0  ;...
		      		        20   -eps       40      -eps      -eps      100    20    0.14   0.01  ];
			
			self.kmax = 4;
			self.p = model.params;
			self.p.phin = 1e-5;
			% Limits based on stability_3d
			% Somewhat tighter ranges
			% EMG ranges based on spectfit output generated by sacha, speciically spec10000113.fit
		end

		function p = p_from_params(self,fitted_params) % Map parameters to point
			p = p_from_params@bt.model.full(self,fitted_params);
			p.beta(:) = fitted_params(7)*fitted_params(6);
		end
		
		function params = params_from_p(self,p) % Map point to parameters
			params = [p.gabcd p.alpha(1) p.beta(1)/p.alpha(1) p.t0 p.emg_a];
		end

		function valid = validate_params(self,pars)
			valid = ~(pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || any(pars > self.limits(2,:) | pars < self.limits(1,:) ));
		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			% DON'T FORGET THAT bt.MODEL.FULL.OBJECTIVE HAS A CHECK FOR BETA/ALPHA<20
			pars(7) = pars(6)*pars(7); % Calculate beta
			[chisq,P] = objective@bt.model.full(self,pars);
		end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[initial_values,prior_pp] = initialize_fit@bt.model.full(self,target_f,target_P);
			initial_values(7) = initial_values(7)/initial_values(6);
			prior_pp = self.uniform_priors();
		end
		
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*(params(:,6).*params(:,7))./((params(:,6)+(params(:,6).*params(:,7))).^2);
		end

	end
end


