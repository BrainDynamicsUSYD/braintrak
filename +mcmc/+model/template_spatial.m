classdef template_spatial < mcmc.model.template
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		p
		n_modes = 6
		n_truncate = 6 % Seems like a useful starting point for cosine variations
		normalization_target
		
		output_x
		output_y
		weights
		emg % This is the actual EMG component as a function of frequency
		emg_scaling
	end

	methods(Abstract)
		set_params(p,pars) % Takes in parameters, modifies the internal p struct appropriately
		% Derived classes must also implement p_from_params
	end

	methods

		function set_electrodes(self,data)
			set_electrodes@mcmc.model.template(self,data)

			% Also calculate the electrode positions
			[self.output_x,self.output_y] = mcmc.electrode_positions(self.electrodes,self.p.Lx,self.p.Lx);
		end

		function [p,p_electrode] = p_from_params(self,fitted_params) 
			% Use the set_params function
			self.set_params(fitted_params);
			p = self.p.copy();
			p.disable_set = false;

			if nargout > 1
				for j = 1:length(self.electrodes)
					p_electrode(j) = p.params_at_xy(self.output_x(j),self.output_y(j));
				end
			end
		end

		
		function params = params_from_p(self,p) % Map point to parameters
			error('Cannot retrieve spatial variations from parameter struct in general')
		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			P = [];
			chisq = NaN;

			% Load params 
			self.set_params(pars);

			% Check stability using is_stable method
			% Will intelligently choose XYZ or Gabcd or Gab stability
			if ~self.p.is_stable()
				return
			end

			% Get nonuniform spectrum
			[~,~,~,~,~,~,~,P] = self.p.spatial_spectrum(self.n_modes,self.n_truncate,self.output_x,self.output_y,self.target_f,false);

			P = shiftdim(P,2);
			P = bsxfun(@plus,P,1e-9*self.p.emg_a*self.emg*2*pi); % Add the EMG component
			P = bsxfun(@rdivide,P,mex_trapz(self.target_f(self.weights>0),P(self.weights>0,:)));
			P = bsxfun(@times,self.normalization_target,P);
			
			sqdiff = (abs(P-self.target_P)./self.target_P).^2; % This is the squared fractional difference
			sqdiff = bsxfun(@times,sqdiff,self.weights);
			chisq = sum(sqdiff(:))/size(P,2); % Failing to divide by the number of electrodes leads to the probability becoming 0 due to numerical precision
		end
		
		function set_cache(self,initial_values)
			% Almost nothing needs to be done, apart from setting the normalization targets
			spec_w = self.target_f(:)*2*pi;
			emg_f = 40;
			self.emg = (spec_w(:)/(2*pi*emg_f)).^2./(1+(spec_w(:)/(2*pi*emg_f)).^2).^2;

			self.weights = self.get_weights(self.target_f);

			if ~isempty(self.target_P)
				% If this cache file is being generated for a specific fit instance
				for j = 1:size(self.target_P,2) % For each electrode
					self.normalization_target(j) = mex_trapz(self.target_f(self.weights>0),squeeze(self.target_P(self.weights>0,j)));
				end
			else
				self.normalization_target = 1;
			end
		end
		
		% TODO: This should probably be removed
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
		end
		
		function prepare_for_fit(self,varargin)
			% For this family of models, skip_fit must be provided manually or else 
			% fitting all parameters is assumed 
			if nargin < 6 || isempty(varargin{5})
				skip_fit = zeros(1,self.n_params);
			else
				assert(length(varargin{5}) == self.n_params);
				self.skip_fit = varargin{5};
			end
				
			prepare_for_fit@mcmc.model.template(self,varargin{1:4});
		end

		% TODO: Placeholders for other types of variations?
		function valid = check_cosine_validity(self,pars,par_idx,amp_idx)
			valid = (pars(par_idx)+abs(pars(amp_idx)))<=self.limits(2,par_idx) && (pars(par_idx)-abs(pars(amp_idx)))>=self.limits(1,par_idx);
		end

	end
end



