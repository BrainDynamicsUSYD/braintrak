classdef nus_mass < bt.model.template
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit

	properties
		kmax
		k0
		k2u
		k2_volconduct

		stab_gamma_prefactor
		spec_gamma_prefactor

		stab_w
		spec_w
		
		p
		p2
		emg
		normalization_target
		stab_Mtot
		spec_Mtot
		compute_Mtot
	end


	methods
		function self = nus_mass(p) % Constructor
			self.name = 'nus_mass';
			self.n_params = 8;
			self.param_names = {'nu_ee','nu_ei','nu_es','nu_se','nu_sr','nu_sn','nu_re','nu_rs'};
			self.param_symbols = {'nu_{ee}','nu_{ei}','nu_{es}','nu_{se}','nu_{sr}','nu_{sn}','nu_{re}','nu_{rs}'};
			self.param_units = {'mV s','mV s','mV s','mV s','mV s','mV s','mV s','mV s'};

			self.n_fitted = 8;
			self.skip_fit = zeros(1,self.n_fitted);

			self.initial_step_size = 1e-4*[1 1 1 1 1 1 1 1];
			self.limits = 1e-3*[ 0.3    -60    0.02  0.03  -32     0.01  0.07  0.02; ...
					60   -0.04   40    48    -0.02   40    6    7];
			
			self.kmax = 4;
			self.k0 = 10; % Volume conduction parameter
			self.p = p.copy();
			self.p.phin = 1e-5;

		end

		function p = p_from_params(self,fitted_params) % Map parameters to point
			p = model.params;
			p.nus = fitted_params;
		end
		
		function params = params_from_p(self,p) % Map point to parameters
			params = [p.nus];
		end

		function valid = validate_params(self,pars)
			valid = true;

			valid = ~any(pars >= self.limits(2,:) | pars <= self.limits(1,:) );

			if ~valid
				return
			end

			self.p2 = self.p.copy();
			self.p2.nus = pars;

			if isempty(self.p2.phia)
				valid = false;
				return
            end

			philim = [ 0.1  0.1  0.1;...
			        60  60  60];

			valid = ~any(self.p2.phia >= philim(2,:) | self.p2.phia <= philim(1,:) );

		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			P = [];
			chisq = NaN;

			self.p2 = self.p.copy();
			self.p2.nus = pars;
             
            try

            catch
                keyboard
            end

		    % stab_w=2*pi*linspace(0,200,5000); % High res array, 0 to 200Hz with 10000 points
		    % [Jee,Jei,Jese,Jesre,Jsrs,Jesn] = get_Jabcd(p,stab_w);
		    % d=((1-1i*stab_w/p.gammae).^2.*(1-Jei)-Jee).*(1-Jsrs)-(Jese + Jesre);      
		    % stab=d(1)>0 && ~any(real(d(2:end))<0 & imag(d(2:end)).*imag(d(1:end-1))<0);

		    % if ~stab
		    % 	return
		    % end

			f = self.target_f(:); % It is critical to start at w=0 for the convolution - but no convolutions anymore
	        w = f*2*pi;
	        dw = w(2)-w(1);
		    [Jee,Jei,Jese,Jesre,Jsrs,Jesn] = get_Jabcd(self.p2,w);
		        
		    P = zeros(size(w));
			%Dew = (1-1i*w./self.p2.gammae).^2;

		    [De,Di,Dr,Ds,Dn] = deal(1);

            %De = Dew + k*p.re.^2;
            %De = Dew;

            A = De.*(Di-Jei)./Di - Jee;

            this_phie = Jesn.*self.p2.phin./(A.*(1-Jsrs) - Jese - Jesre);

            P = abs(this_phie).^2; 
          
	        % P = P.*dk.^2; % Multiply by dk then sum to get the integral over k
	        P = P(:).'*2*pi; % Convert to P(f)
    	
		    P = P./utils.mex_trapz(self.target_f(self.weights>0),P(self.weights>0));
            P = self.normalization_target*P(:);

		    sqdiff = (abs(P-self.target_P)./self.target_P).^2; % This is the squared fractional difference
		    chisq = sum(sqdiff(:).*self.weights(:));
		end
		
		function prepare_for_fit(self,varargin)
			% argument list: target_f,target_P,initial_values,prior_pp,skip_fit
			if nargin < 6 || isempty(varargin{5})
				self.skip_fit(:) = false;
			else
				assert(length(varargin{5}) == self.n_params);
				self.skip_fit = varargin{5};
			end
			
			prepare_for_fit@bt.model.template(self,varargin{1:4});
		end			
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			initial_values =  self.p.nus;
			prior_pp = self.uniform_priors();
		end
	
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz = randn(size(params,1),3);
		end
		
		function set_cache(self,initial_values) % Calculate the cache - takes in skip_fit
			self.weights = self.get_weights(self.target_f);
			if ~isempty(self.target_P)
				% If this cache file is being generated for a specific fit instance
				self.normalization_target = utils.mex_trapz(self.target_f(self.weights>0),self.target_P(self.weights>0));
			end
		end

	end
end

function [Jee,Jei,Jese,Jesre,Jsrs,Jesn] = get_Jabcd(p,w)
    L = @(idx) 1./((1-1i*w/p.alpha(idx)).*(1-1i*w/p.beta(idx)));
    M = exp(1i*w*p.taues); % Assume Mes = Mse

    % It is important to use the G_ab if they are available
    % because the G_ab define a value for G_sn which allows the
    % normalization with NeuroField to be correct.
    % The 5 gains do not specify the normalization of the power spectrum
    if isempty(p.gab)
        Jee = L(1).*p.gabcd(1);
        Jei = L(2).*p.gabcd(2);
        Jese = L(3).*L(4).*p.gabcd(3).*M.*M;
        Jesre = L(3).*L(5).*L(7).*p.gabcd(4).*M.*M;
        Jsrs = L(5).*L(8).*p.gabcd(5);
        Jesn = L(3).*L(6); 
    else
        % Use the 8 gains
        Jee = L(1).*p.gab(1);
        Jei = L(2).*p.gab(2);
        Jes = L(3).*p.gab(3).*M;
        Jse = L(4).*p.gab(4).*M;
        Jsr = L(5).*p.gab(5);
        Jsn = L(6).*p.gab(6);
        Jre = L(7).*p.gab(7).*M;
        Jrs = L(8).*p.gab(8);

        Jese = Jes.*Jse;
        Jesre = Jes.*Jsr.*Jre;
        Jsrs = Jsr.*Jrs;
        Jesn = Jes.*Jsn;
    end
end
