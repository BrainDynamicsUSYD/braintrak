classdef dummy < mcmc.model.template

	properties
		means 
		std 
	end


	methods



		function self = dummy() % Constructor
			self.name = 'dummy';
			self.electrodes = {'Electrode'};
			self.n_params = 4;
			self.param_names = {'Unimodal','Unimodal','Exponential','Bimodal'};
			self.param_symbols =  {'Unimodal','Unimodal','Exponential','Bimodal'};
			self.param_units = {'','','',''};

			self.n_fitted = 4;
			self.skip_fit = zeros(1,self.n_fitted);

			self.initial_step_size = 0.1*[1 1 1 1];
		   	self.limits = [ -20 -20 -20 -20  ;...
		      		        20  20  20  20  ];
			
	
		end

		function set_electrodes(self,a)
		end

		function [initial_params,initial_pp] = initialize_fit(self,a,b)
			initial_params = ones(1,4);
			initial_pp = self.uniform_priors();
		end

		% function prior_pp = uniform_priors(self)
		%     for j = 1:size(self.limits,2)
		%         a=linspace(self.limits(1,j),self.limits(2,j),200);
		%         b=ones(200,1);
		%         prior_pp.x(:,j) = a(:);
		%         prior_pp.y(:,j) = b(:);
		%         prior_pp.ndx(j) = 1/(a(2)-a(1));
		%     end
		% end

		function prepare_for_fit(self,a,b,c,d,e);
        end

		function p = p_from_params(self,pars)
			p = model.params;
		end

		function params = params_from_p(self)
			params = [];
		end

		function validate(self,a)
		end

		function p = eval_priors(self,xi)
			p = 1;
		end

		function objective(self,x);
		end

		%[fit_data.fitted_chisq,~,likelihood,fit_data.fitted_P]

		function [chisq,prob,likelihood,fitted_P] = probability(self,pars)
			chisq = 1;
			prob = NaN;
			if any(pars < self.limits(1,:) | pars > self.limits(2,:))
				chisq = NaN;
				return;
			end

			means = [1,-3,2,-5 20];
			stds = [4, 12, NaN,2, 3];

			p = zeros(1,4);

			p(1) = self.g(pars(1),means(1),stds(1));
			p(2) = self.g(pars(2),means(2),stds(2));
			lambda = means(3).^-1;
			p(3) = lambda*exp(lambda*pars(3)).*(pars(3)>0);
			p(4) = self.g(pars(4),means(4),stds(4)) + self.g(pars(4),means(5),stds(5));

			likelihood =  prod(p); % Probability of the state
			prob = likelihood*self.eval_priors(pars);
			fitted_P = 1;
		end

		function xyz = get_xyz(self,pars)
			xyz = repmat([0.1 0.1 0.1],size(pars,1),1);
		end


		function v = g(self,x,m,s)
			v = exp(-(x-m)^2/(2*s^2));
		end

	end
end


% function [chisq,prob,likelihood,P] = probability(self,pars)
% 	P = [];
% 	chisq = NaN;
% 	likelihood = NaN;
% 	prob = NaN;

% 	% Use the base class to validate params
% 	% this will occur before any adapter code
% 	if ~self.validate_params(pars)
% 		return
% 	end

% 	[chisq,P] = self.objective(pars);
% 	if isfinite(chisq)
% 		likelihood = exp(-chisq/2);
% 		prob = likelihood*self.eval_priors(pars);
% 	else
% 		likelihood = 0;	
% 		prob = 0;
% 	end
% end