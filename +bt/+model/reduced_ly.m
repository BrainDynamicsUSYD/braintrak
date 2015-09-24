classdef reduced_ly < bt.model.template
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit
	% Note that this does sum over k

	properties
		kmax
		k0
		k2u
		k2_volconduct

		stab_w
		spec_w
		p
		stab_gamma_prefactor
		spec_gamma_prefactor
		emg
		normalization_target
		stab_Mtot
		spec_Mtot
		compute_Mtot
	end


	methods
		function self = reduced_ly() % Constructor
			self.name = 'model_reduced_ly';
			self.n_params = 7;
			self.param_names = {'X','Y','Z','Alpha','Beta','t0','EMGa'};
			self.param_symbols = {'X','Y','Z','\alpha','\beta','t_0','A_{EMG}'};
			self.param_units = {'','','','s^{-1}','s^{-1}','ms',''};
			
			self.n_fitted = 7;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.05.*[0.2 0.2 0.2 5 10 0.02 0.05];
		   	self.limits = [ eps  -1   eps       10      100    0.075   0  ;...
		      		        1     1     1      100      800    0.14    5  ];

		    self.kmax = 4; % Number of modes to simulate
			self.p = model.params;
			self.p.phin = 1e-5;
			% Limits based on stability_3d and for reduced model, classic tent
			% Somewhat tighter ranges
			% EMG ranges based on spectfit output generated by sacha, speciically spec10000113.fit
		end

		function p =  p_from_params(self,fitted_params) % Map parameters to point
			p = model.params;
			p.xyz = fitted_params(1:3);
			p.alpha(:) = fitted_params(4);
			p.beta(:) = fitted_params(5);
			p.t0 = fitted_params(6);
			p.taues = p.t0/2;
			p.tause = p.t0/2;
			p.emg_a = fitted_params(7);
		end
		
		function params = params_from_p(self,p) % Map point to parameters
			params = [p.xyz p.alpha(1) p.beta(1) p.t0 p.emg_a];
		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			chisq = NaN;
			P = [];
			if any(pars > self.limits(2,:) | pars < self.limits(1,:))
				return
			end

		    if self.compute_Mtot % If this flag has been set to false, then these fields are already present
		        % Should be fine to let that just throw an undefined variable error later if this isn't the case
		       self.stab_Mtot = exp(1i*self.stab_w*pars(6));
		       self.spec_Mtot = exp(1i*self.spec_w*pars(6));
		    end

			stab_L = 1./((1-1i*self.stab_w/pars(4)).*(1-1i*self.stab_w/pars(5)));
		    spec_L = 1./((1-1i*self.spec_w/pars(4)).*(1-1i*self.spec_w/pars(5)));
		    zprime = (pars(4)+pars(5))^2/(pars(4)*pars(5))*pars(3); % Effectively Gsrs

		    % cf. James corticothalamic dynamics eq 9
		    d = ((self.stab_gamma_prefactor) - pars(1)).*(1+stab_L.*stab_L.*zprime) - stab_L.*stab_L.*pars(2).*(1+zprime).*self.stab_Mtot;
		    stab=d(1)>0 && ~any(real(d(2:end))<0 & imag(d(2:end)).*imag(d(1:end-1))<0);

		    if ~stab
		    	return
		    end

		    % And the spectrum. At k = 0
		    one_plus_zprime = (1+zprime.*spec_L.*spec_L);

		    P = zeros(size(self.spec_w));
		    q2re2 = (self.spec_gamma_prefactor - pars(1) -(spec_L.*spec_L.*pars(2)*(1+zprime))./one_plus_zprime.*self.spec_Mtot);
		    k2u = self.k2u;
		    k2_volconduct = self.k2_volconduct;
		    re2 = self.p.re.^2;

		    for j = 1:size(k2u,1)
		        P = P + k2u(j,2).*abs((spec_L.*spec_L)./(one_plus_zprime.*(k2u(j,1)*re2+q2re2))).^2 * k2_volconduct(j); % For use with the efficient way
		    end

		    P = P + 1e-1*pars(7)*self.emg; % Add the EMG component
		    
		    P = P./mex_trapz(self.target_f(self.weights>0),P(self.weights>0));
		    P = self.normalization_target*P;
		    sqdiff = (abs(P-self.target_P)./self.target_P).^2; % This is the squared fractional difference
		    chisq = sum(sqdiff(:).*self.weights(:));
		end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data] = db_fit.quick_fit(self,target_f,target_P);
			a =  db_data.xyz(idx,:);
			p = model.params(db_data.iswake(idx));

			initial_values =  [a(1) a(2) a(3) p.alpha(1) p.beta(1) p.t0 0];

			prior_pp = self.uniform_priors();
		end
		
		function set_cache(self,initial_values)
			% Get the k variables
			Lx = 0.5; % linear dimensions of cortex (assume square)
			self.k0 = 10; % Volume conduction parameter
			dk = 2*pi/Lx;
			m_rows = -self.kmax:self.kmax; 
			n_cols = -self.kmax:self.kmax; 
			[kxa,kya] = meshgrid(dk*m_rows,dk*n_cols); 
			k2 = kxa.^2+kya.^2;
			k2u = unique(k2(:));
			self.k2u = [k2u histc(k2(:),k2u)]; % First col is the k2 value, second is how many times it appeared
			self.k2_volconduct = exp(-k2u(:,1)/self.k0^2);


			% Cache for the mcmc fit code
			% Lower resolution, dynamic alpha and beta
			self.stab_w=2*pi*linspace(0,100,5000); % High res array, 0 to 100Hz with 10000 points
			self.spec_w = self.target_f*2*pi;

			self.stab_w = self.stab_w(:);
			self.spec_w = self.spec_w(:);

			self.stab_gamma_prefactor = (1-1i*self.stab_w/self.p.gammae).^2;
			self.spec_gamma_prefactor = (1-1i*self.spec_w/self.p.gammae).^2;

			self.weights = self.get_weights(self.target_f);


			emg_f = 40;
			self.emg = (self.spec_w/(2*pi*emg_f)).^2./(1+(self.spec_w/(2*pi*emg_f)).^2).^2;
			
			if ~isempty(self.target_P)
				% If this cache file is being generated for a specific fit instance
				self.normalization_target = mex_trapz(self.target_f(self.weights>0),self.target_P(self.weights>0));
			end


			if self.skip_fit(6) == 2 && ~isempty(initial_values) % If t0 is not being fitted, then the value is of course fixed
				self.stab_Mtot = exp(1i*self.stab_w*initial_values(6));
				self.spec_Mtot = exp(1i*self.spec_w*initial_values(6));
				self.compute_Mtot = false;
			else
				self.compute_Mtot = true; 
			end
		end
		
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1);
			xyz(:,2) = params(:,2);
			xyz(:,3) = params(:,3);
		end
		
		function prepare_for_fit(self,varargin)
			% argument list: target_f,target_P,initial_values,prior_pp,skip_fit
			if nargin < 6 || isempty(varargin{5})
				self.default_t0_skip_fit(varargin{1},varargin{2},6);
			else
				assert(length(varargin{5}) == self.n_params);
				self.skip_fit = varargin{5};
			end
			
			prepare_for_fit@bt.model.template(self,varargin{1:4});
		end

	end
end