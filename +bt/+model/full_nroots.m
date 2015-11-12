classdef full_nroots < bt.model.template
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit

	properties
		p
		ve_test_values
	end


	methods
		function self = full_nroots() % Constructor
			self.name = 'model_full_nroots';
			self.n_params = 8;
			self.param_names = {'nuee','nuei','nues','nuse','nusr','nusn','nure','nurs'};
			self.param_symbols = {'nu_ee','nu_ei','nu_es','nu_se','nu_sr','nu_sn','nu_re','nu_rs'};
			self.param_units = {'mV\,s','mV\,s','mV\,s','mV\,s','mV\,s','mV\,s','mV\,s','mV\,s'};

			self.n_fitted = 8;
			self.skip_fit = zeros(1,self.n_fitted);

			self.initial_step_size = 5e-3*[10 1 1 1 1 1 1 1];
		   	self.limits = 1e-3*[ 0  -60        0     0      -32      0      0   0  ;...
		      		            340    0       40       48       0      40    64    64  ];
			
			self.p = model.params;
			self.ve_test_values = -0.04:0.00002:0.04;

		end

		function p = p_from_params(self,fitted_params) % Map parameters to point
			p = model.params;
			p.nus = fitted_params;
		end
		
		function params = params_from_p(self,p) % Map point to parameters
			params = p.nus;
		end

		function valid = validate_params(self,pars)
			valid = ~(any(pars >= self.limits(2,:) | pars <= self.limits(1,:) ));
		end


		function objective(self,x);
		end

		function [chisq,prob,likelihood,fitted_P] = probability(self,pars) % Calculate the objective
			P = [];
			chisq = NaN;
			likelihood = 0;
			fitted_P = 1;
			prob = 0;

			if ~self.validate_params(pars)
				return
			end

			%[ves,discard] = utils.allroots(@model.mex_ve_root,self.ve_test_values,pars,self.p.theta,self.p.sigma,self.p.qmax);
			[ves,discard] = utils.allroots(model.mex_ve_root(self.ve_test_values,pars,self.p.theta,self.p.sigma,self.p.qmax),self.ve_test_values);

			% chisq should be close to zero (but not quite) 
			% when there are 5 roots
			chisq = 0.1+5-length(ves); % The number of roots is chisq - 5.1
			likelihood =  exp(-chisq*2); % Probability of the state
			prob = likelihood;
		end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			pts = load('corticothalamic-model/example_parameters.mat');
			initial_values = pts.ec.nus;
			prior_pp = self.uniform_priors();
		end
		
		function set_cache(self,initial_values) % Calculate the cache - takes in skip_fit
			return
		end
		
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz = [NaN NaN NaN];
		end
		
		function prepare_for_fit(self,varargin);
        end

	end
end


