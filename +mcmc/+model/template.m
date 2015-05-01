classdef (Abstract) template < matlab.mixin.Copyable
    
	% These properties are SetAccess protected to ensure that the MCMC functions
	% are properly independent of the choice of model
	properties(SetAccess = protected) % Public access, but can only be edited via the interface
		name
		n_params
		param_names
		param_symbols
		param_units

		n_fitted % Expect to change during fitting. Number of parameters being fitted. this could be different if skip_fit changes
		initial_step_size % Standard deviations for initial chain (before adaptive covariance matrix is calculated)
		limits % Upper and lower bounds on parameter values
		skip_fit % Expect to change during fitting. 1 means distribution is fixed, 2 means parameter AND distribution are fixed
		electrodes % So that the model can choose its behaviour based on the electrodes
	
		target_f % Frequencies present in the data 
		target_P % Experimental power spectrum
		prior_pp % Shape of the prior distribution (a struct array, one for each parameter)

	end

   	methods(Abstract)
		p_from_params(self,params) % Map parameters to point
		params_from_p(self,p) % Map point to parameters
		objective(self,pars) % Calculate the objective and the spectrum. Note that this syntax forces prepare_for_fit to be run
		initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
		get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
	end


	methods

		function [target_f,P,stab] = spectrum(self,pars,target_f)

			if nargin < 3
				target_f = linspace(1,45,100);
			end
			self.target_f = target_f(:);
			self.target_P = ones(size(self.target_f))./(target_f(2)-target_f(1));
			self.set_cache(pars);
			[chisq,P] = self.objective(pars);
			stab = isfinite(chisq);
			P = P./utils.mex_trapz(target_f,P);
		end

		function default_t0_skip_fit(self,target_f,target_P,t0_index)
			% Default skip fit behaviour for t0
			% Replace this function internally if different behaviour is required

			self.skip_fit = zeros(1,self.n_params);
			assert(isvector(target_P)); % This function should not be used with multiple electrodes
			temp_analysis = model.get_spec_analysis(target_f(:),target_P(:));
			if temp_analysis.alpha_str < 5 % Only fit t0 if alpha peak is present
				self.skip_fit(t0_index) = 1; % Change this to false to enable fitting t0 at all times
			else
				self.skip_fit(t0_index) = 0;
			end
		end

		function prepare_for_fit(self,target_f,target_P,initial_values,prior_pp)
			self.n_fitted = sum(~self.skip_fit);
			self.target_f = target_f(:);
			self.target_P = target_P;
			self.prior_pp = prior_pp;
			self.set_cache(initial_values);
		end

		function validate(self)
			% This function is the first thing called by chain()
			% and provides a last opportunity to detect any problems with the
			% internal state of the model before model.probability is called
			% Perform any checks/assertions required immediately before fitting
			if any(self.target_f==0)
				error('Weighting will only work with nonzero frequencies');
			end
		end

		function w = get_weights(self,target_f) % Default weighting function
			w = ones(size(target_f));
			w = target_f.^-1; % Decades equally weighted
			%w(target_f < 1) = 0; % Ignore frequencies below 1Hz 
			w = w(:);
		end

		function set_electrodes(self,data)
			if isstruct(data)
				self.electrodes = data.colheaders;
			elseif iscell(data)
				self.electrodes = data;
			else
				self.electrodes = {data};
			end
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

		% This is called by probability
		% Overload this to perform more sophisticated checks
		function valid = validate_params(self,pars)
			valid = ~any(pars > self.limits(2,:) | pars < self.limits(1,:));
		end

		function [chisq,prob,likelihood,P] = probability(self,pars)
			P = [];
			chisq = NaN;
			likelihood = NaN;
			prob = NaN;

			% Use the base class to validate params
			% this will occur before any adapter code
			if ~self.validate_params(pars)
				return
			end

			[chisq,P] = self.objective(pars);
			if isfinite(chisq)
				likelihood = exp(-chisq/2);
				prob = likelihood*self.eval_priors(pars);
			else
				likelihood = 0;	
				prob = 0;
			end
		end

		function p = eval_priors(self,xi)
			% Evaluates all priors and returns the product of them
			% Note that the proposed values must have been checked to lie within the allowed limits
			% BEFORE this function is called
			xi = xi - self.prior_pp.x(1,:);
			fxi = floor(xi.*self.prior_pp.ndx)+1;
			Yi = (fxi-xi.*self.prior_pp.ndx).*diag(self.prior_pp.y(fxi,:))'+(1-fxi+xi.*self.prior_pp.ndx).*diag(self.prior_pp.y(fxi+1,:))';
			p = prod(Yi);
		end

		function posterior_pp = make_posterior(self,out,lim)
			% Given the output run (npts x 8 params)
			% Return the posterior distribution in pp format
			if nargin < 3 || isempty(lim)
				lim = self.limits;
			end

			for j = 1:size(lim,2)
				%x1 = linspace(lim(1,j),lim(2,j),15); % Histogram bins
				%posterior_pp.x(:,j) = linspace(lim(1,j),lim(2,j),50); % Final resolution
				posterior_pp.x(:,j) = linspace(lim(1,j),lim(2,j),50);
				[posterior_pp.y(:,j)] = ksdensity(out(:,j),posterior_pp.x(:,j));

				%y1 = hist(out(:,j),x1);
				%posterior_pp.y(:,j) = csaps(x1,y1,1,posterior_pp.x(:,j)) 
				%posterior_pp.y(:,j) = smooth(posterior_pp.y(:,j),0.3,'rloess');
				posterior_pp.y(:,j) = posterior_pp.y(:,j)/trapz(posterior_pp.x(:,j),posterior_pp.y(:,j));
				posterior_pp.ndx(j) = 1/(posterior_pp.x(2,j)-posterior_pp.x(1,j));
			end
		end

		function post = xyz_posterior(self,out)
			xyz_vals = self.get_xyz(out);
			xyz_lim = [0  -1     0 ; 1   1   1.3];
			post = self.make_posterior(xyz_vals,xyz_lim);
		end

		function present = xyz_present(self)
			% Returns [true true true] if there is a parameter
			% with name 'X','Y', or 'Z' respectively
			% e.g. [1 1 1] for reduced model, [0 1 1] for reduced_l
			% and [0 0 0] for full model
			present = [any(strcmp(self.param_names,'X')) any(strcmp(self.param_names,'Y')) any(strcmp(self.param_names,'Z'))];
		end
	end
end

