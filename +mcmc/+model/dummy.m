classdef dummy < handle
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit

	properties
		name 
		electrodes 
		n_params 
		param_names 
		means 
		std 
		n_fitted 
		skip_fit 
		initial_step_size 
		limits 
	end


	methods



		function self = dummy() % Constructor
			self.name = 'dummy';
			self.electrodes = 1;
			self.n_params = 5;
			self.param_names = {'A','B','C','D','E'};

			self.means = [5 4 3 2 1];

			self.std = [1 2 3 4 5];

			self.n_fitted = 5;
			self.skip_fit = zeros(1,self.n_fitted);

			self.initial_step_size = 0.2*[1 1 1 1 1];
		   	self.limits = [ -20 -20 -20 -20 -20 ;...
		      		        20  20  20  20  20  ];
			
	
		end

		function posterior_pp = make_posterior(self,out)
			% Given the output run (npts x 8 params)
			% Return the posterior distribution in pp format
			for j = 1:size(self.limits,2)
				posterior_pp.x(:,j) = linspace(self.limits(1,j),self.limits(2,j),200);
				posterior_pp.y(:,j) = hist(out(:,j),posterior_pp.x(:,j));
				posterior_pp.y(:,j) = smooth(posterior_pp.y(:,j),15);
				posterior_pp.y(:,j) = posterior_pp.y(:,j)/trapz(posterior_pp.x(:,j),posterior_pp.y(:,j));
				posterior_pp.ndx(j) = 1/(posterior_pp.x(2,j)-posterior_pp.x(1,j));
			end
		end

		function set_electrodes(self,a)
		end

		function [initial_params,initial_pp] = initialize_fit(self,a,b)
			initial_params = ones(size(self.means));
			initial_pp = self.uniform_priors();
		end

		function prior_pp = uniform_priors(self)
		    for j = 1:size(self.limits,2)
		        a=linspace(self.limits(1,j),self.limits(2,j),200);
		        b=ones(200,1);
		        prior_pp.x(:,j) = a(:);
		        prior_pp.y(:,j) = b(:);
		        prior_pp.ndx(j) = 1/(a(2)-a(1));
		    end
		end

		function prepare_for_fit(self,a,b,c,d,e);
		end


		function validate(self,a)
		end

		function [chisq,current_posterior] = probability(self,x)
			chisq = 1;
			current_posterior = NaN;
			if any(x < self.limits(1,:) | x > self.limits(2,:))
				chisq = NaN;
				return;
			end
			p = zeros(length(self.means),1);
			for j = 1:length(p)
				p(j) = self.g(x(j),self.means(j),self.std(j));
			end
			current_posterior = prod(p);
		end


		function v = g(self,x,m,s)
			v = exp(-(x-m)^2/(2*s^2));
		end

	end
end


