classdef spatial_t0_2d < mcmc.model.template
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit

	properties
		Lx
		dk
		dx
		grid_x
		grid_y
		grid_kx
		grid_ky
		grid_K2
		K2_volconduct
		K2re2_diag % To add to the LHS after taking the FFT and convolution masking
		T_mask % The 2D circulant convolution matrix index mask

		output_x % The electrode positions to evaluate the spectrum at
		output_y % These are set by the local set_electrodes() method

		stab_w
		spec_w
		p
		stab_gamma_prefactor
		spec_gamma_prefactor
		weights
		k2_volconduct
		emg
		normalization_target
		stab_Mtot
		spec_M
		compute_Mtot
	end


	methods
		function self = spatial_t0_2d(n_modes) % Constructor
			self.name = 'model_spatial_t0_2d';
			self.n_params = 10;
			self.param_names = {'Gee','Gei','Gese','Gesre','Gsrs','Alpha','Beta','t0','EMGa','t0_amp'};


			self.n_fitted = 10;
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

			% Initialize the spatial variables
			p = model.params;
			if nargin < 1 || isempty(n_modes)
				n_modes = 4; % Number of points to use for the spatial variations and for the spatial modes
			end

			self.Lx = 0.5; % linear dimensions of cortex
			[~,self.grid_kx,self.grid_ky,self.grid_x,self.grid_y,~,self.dk,self.dx] = get_frequencies(zeros(n_modes,n_modes,1),1,self.Lx,self.Lx);

			self.grid_kx = ifftshift(self.grid_kx);
			self.grid_ky = ifftshift(self.grid_ky); % Set up these arrays so that the forward transform doesn't need shifting
			self.grid_K2 = self.grid_kx.^2+self.grid_ky.^2;
			k0 = 10; % Volume conduction parameter
			self.K2_volconduct = exp(-self.grid_K2/k0^2);
			self.K2_volconduct = self.K2_volconduct(:)'; % Row vector so bsxfun works properly in objective()

			Qm = fft(eye(n_modes));
			Q = kron(Qm,Qm);
			Lm = fft2(reshape(1:n_modes^2,[n_modes,n_modes])); % Note - need to undo the  applied to ap initially
			self.T_mask = round(Q'*diag(Lm(:))*Q/numel(Lm)); 
			self.K2re2_diag = diag(self.grid_K2(:)*p.re^2);
		end

		function set_electrodes(self,data)
			set_electrodes@mcmc.model.template(self,data)

			% Also calculate the electrode positions
			[self.output_x,self.output_y] = mcmc.electrode_positions(self.electrodes,self.Lx,self.Lx);
		end

		function p = p_from_params(self,fitted_params) % Map parameters to point
			p = model.params;
			p.gabcd = fitted_params(1:5);
			p.alpha(:) = fitted_params(6);
			p.beta(:) = fitted_params(7);
			p.t0 = fitted_params(8);
			p.taues = p.t0/2;
			p.tause = p.t0/2;
			p.spatial_t0 = @(x,y) fitted_params(8)+fitted_params(10)*cos(y*2*pi/self.Lx);
			p.emg_a = fitted_params(9);
		end
		
		function params = params_from_p(self,p) % Map point to parameters
			error('Cannot retrieve spatial variations from parameter struct!')
			%params = [p.gabcd p.alpha(1) p.beta(1) p.t0 p.t0_amp p.emg_a];
		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			P = [];
			chisq = NaN;
			if pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || pars(7)/pars(6) > 20 || any(pars > self.limits(2,:) | pars < self.limits(1,:) )
				return
			end

			%------------------- CHECK STABILITY

			if self.compute_Mtot % If this flag has been set to false, then these fields are already present
				% Calculate the grid of t0 values - specific to this model!
				self.stab_Mtot = exp(1i*self.stab_w*pars(8)); % Use the center value
			end
			

			stab_L = 1./((1-1i*self.stab_w/pars(6)).*(1-1i*self.stab_w/pars(7)));
			d=(self.stab_gamma_prefactor.*(1-stab_L.*pars(2))-stab_L.*pars(1)).*(1-stab_L.*stab_L.*pars(5))-stab_L.*stab_L.*self.stab_Mtot.*(pars(3) + stab_L.*pars(4));      

			stab=d(1)>0 && ~any(real(d(2:end))<0 & imag(d(2:end)).*imag(d(1:end-1))<0);

			if ~stab
				return
			end

			%------------------- CALCULATE SPECTRUM

			if self.compute_Mtot % If this flag has been set to false, then these fields are already present
				% Calculate the grid of t0 values - specific to this model!
				spec_t0 = pars(8)+pars(10)*cos(self.grid_y*2*pi/self.Lx); % Mixing y and x assumes a square cortex
				self.spec_M = exp(1i*bsxfun(@times,self.spec_w,spec_t0/2)); % This is gridded
			end

			% This could also be spatially and frequency dependent at a later date
			% In which case, it would be the same size as self.spec_M
			% All that needs to change is some of the bsxfun calls below will become 
			% simple element-wise operations
			spec_L = 1./((1-1i*self.spec_w/pars(6)).*(1-1i*self.spec_w/pars(7)));
			Jei_oneminus = 1-spec_L.*pars(2);
			Jsrs_oneminus = 1-spec_L.*spec_L.*pars(5);
			Afull = (spec_L.*spec_L.*pars(3) + spec_L.*spec_L.*spec_L.*pars(4))./(Jei_oneminus.*Jsrs_oneminus);
			Afull = bsxfun(@times,Afull,(self.spec_M.*self.spec_M));
			Afull = bsxfun(@minus,self.spec_gamma_prefactor - (spec_L.*pars(1))./Jei_oneminus,Afull);
			Bfull = bsxfun(@times,(spec_L.*spec_L.*1)./((Jei_oneminus).*(Jsrs_oneminus)),self.spec_M); % Assume Gesn = 1

			% Then take the 2D FFT of both of these quantities
			Afull = fft(Afull,[],1);
			Afull = fft(Afull,[],2);
			Afull = Afull/numel(self.grid_K2); % verified correct normalization
			Bfull = fft(Bfull,[],1);
			Bfull = fft(Bfull,[],2);
			Bfull = Bfull/numel(self.grid_K2); % verified correct normalization

			% Preallocate the output matrices
			P = zeros([length(self.spec_w) length(self.output_x)]);

			% -- LOOP OVER FREQUENCY
			for j = 1:size(Afull,3)
				A_temp = Afull(:,:,j);
				B_temp = Bfull(:,:,j);
				A = A_temp(self.T_mask);
				A = A+self.K2re2_diag;
				B = B_temp(self.T_mask);
				M = A\B; % inv(A)*B == A\B

				mm = M*M';
				mm = bsxfun(@times,mm,self.K2_volconduct);

				if j == 1 || j == 2
					pt = zeros(size(self.output_x));
					mm_filter = abs(mm) > 0;
					kgx = bsxfun(@minus,self.grid_kx(:),self.grid_kx(:)');
					kgy = bsxfun(@minus,self.grid_ky(:),self.grid_ky(:)');
					kgx1 = kgx(mm_filter);
					kgy1 = kgy(mm_filter);
					exponent_cache = zeros(length(kgx1),numel(self.output_x));
					for j_1 = 1:numel(self.output_x)
						exponent_cache(:,j_1) = exp(1i*kgx1*self.output_x(j_1) + 1i*kgy1*self.output_y(j_1));
					end
				end

				mm_terms = mm(mm_filter);

				for j_1 = 1:numel(self.output_x) % For each position
					tmp = exponent_cache(:,j_1) .* mm_terms;
					P(j,j_1) = abs(sum(tmp(:)));
				end

			end
			
			% This normalization block is for comparing with harness.m and comparing with NeuroField but
			% can be skipped when actually fitting as they get rolled into the experimental normalization
			p.phin = 1.0000e-05; % Copied from wake_latest while testing
			phin = p.phin*ones(numel(self.grid_K2),1);
			phin_multiply = sum(abs(phin).^2)*self.dk*self.dk; % Get phin(w) from phin(k,w)
			P = P.*phin_multiply/numel(self.grid_kx)*self.dk*self.dk/2/pi/2/pi; 
			P = P.*2*pi;

			% Otherwise, add EMG and continue
			P = bsxfun(@plus,P,1e-9*pars(9)*self.emg*2*pi); % Add the EMG component
			P = bsxfun(@rdivide,P,mex_trapz(self.target_f(self.weights>0),P(self.weights>0,:)));
			P = bsxfun(@times,self.normalization_target,P);

			sqdiff = (abs(P-self.target_P)./self.target_P).^2; % This is the squared fractional difference
			sqdiff = bsxfun(@times,sqdiff,self.weights);
			chisq = sum(sqdiff(:));
		end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,squeeze(target_P(:,1))); % Fit the first (or only) spectrum given
			% loglog(target_f,squeeze(target_P(:,4)))
			% hold on
			% loglog(f1,P1)
			% keyboard
			a =  db_data.gab(idx,:);
			p = model.params(db_data.iswake(idx));

			initial_values =  [a(1) a(2) a(3)*a(4) a(3)*a(5)*a(7) a(5)*a(8) p.alpha(1) p.beta(1) p.t0 0 0];

			prior_pp = self.uniform_priors();
		end
		
		function set_cache(self,initial_values)
			% Set up standard frequency variables
			self.stab_w=2*pi*linspace(0,100,5000); % High res array, 0 to 100Hz with 10000 points
			self.spec_w(1,1,:) = self.target_f*2*pi;
			self.stab_w = self.stab_w(:);
			self.p = model.params;
			self.stab_gamma_prefactor = (1-1i*self.stab_w/self.p.gammae).^2;
			self.spec_gamma_prefactor = (1-1i*self.spec_w/self.p.gammae).^2;
			self.weights = self.get_weights(self.target_f);
			emg_f = 40;
			self.emg = (self.spec_w(:)/(2*pi*emg_f)).^2./(1+(self.spec_w(:)/(2*pi*emg_f)).^2).^2;
			
			% Calculate the normalization target for each electrode
			if ~isempty(self.target_P)
				% If this cache file is being generated for a specific fit instance
				for j = 1:size(self.target_P,2) % For each electrode
					self.normalization_target(j) = mex_trapz(self.target_f(self.weights>0),squeeze(self.target_P(self.weights>0,j)));
				end
			else
				self.normalization_target = 1;
			end

			if self.skip_fit(8) == 2 && ~isempty(initial_values) % If t0 is not being fitted, then the value is of course fixed
				% Calculate the grid of t0 values - specific to this model!
				self.stab_Mtot = exp(1i*self.stab_w*initial_values(8)); % Use the center value
				spec_t0 = initial_values(8)+initial_values(10)*cos(self.grid_y*2*pi/self.Lx);
				self.spec_M = exp(1i*bsxfun(@times,self.spec_w,spec_t0/2));
				self.compute_Mtot = false;
			else
				self.compute_Mtot = true; 
			end

		end
		
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
		end
		
		function prepare_for_fit(self,varargin) % Use target_f and target_P to decide if skip_fit should be changed
			% prepare_for_fit needs a friendly interface for users to quickly get the spectrum
			% In general, all arguments are provided by fit()
			if isempty(self.output_x)
				error('Electrode positions were not set correctly. Must run set_electrodes() earlier');
			end

			target_f = varargin{1};
			target_P = varargin{2};

			if nargin < 6 || isempty(varargin{5})
				skip_fit = zeros(1,10);
				alphastr = zeros(1,size(target_P,3));
				
				for j = 1:size(target_P,2)
					temp_analysis = get_spec_analysis(target_f(:),squeeze(target_P(:,j)));
					alphastr(j) = temp_analysis.alpha_str;
					%fprintf('Alpha strength: %.2f\n',temp_analysis.alpha_str);
				end
				
				% Automatically decide fitting
				if sum(alphastr < 5)/size(target_P,2) > 0.5 % If more than half the electrodes have no alpha
					skip_fit([8,10]) = 2; % Change this to false to enable fitting t0 at all times
				else
					skip_fit([8,10]) = 0;
				end
				
			else
				assert(length(varargin{5}) == self.n_params);
				self.skip_fit = varargin{5};
			end
				
			prepare_for_fit@mcmc.model.template(self,varargin{1:4});
		end

	end
end



