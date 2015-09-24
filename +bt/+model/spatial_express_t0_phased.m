classdef spatial_express_t0_phased < bt.model.template_spatial
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
	end

	methods
		function self = spatial_express_t0_phased(n_modes) % Constructor
			self.name = 'spatial_express_t0_phased';
			self.n_params = 11;
			self.param_names = {'Gee','Gei','Gese','Gesre','Gsrs','Alpha','Beta','t0','EMGa','t0_amp','t0_phase'};
			self.param_symbols = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta','t_0','A_{EMG}','t_{0amp}','t_{0\phi}'};
			self.param_units = {'','','','','','s^{-1}','s^{-1}','ms','','ms','m'};

			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.05.*[2 2 2 2 2 5 10 0.02 0.01 0.05 0.2];
			self.limits = [
				eps       	20     	  %Gee
				-20     	-eps      %Gei
				eps       	40     	  %Gese
				-40     	-eps      %Gesre
				-14     	-eps      %Gsrs
				10      	100    	  %Alpha
				100     	800    	  %Beta
				75e-3  		140e-3 	  %t0
				0       	5      	  %EMGa
				-20e-3      20e-3     %t0_amp
				0 			0.5		  % t0_phase
				].';

			self.p = model.params();
			self.p.disable_set = true;
		end

		function w = get_weights(self,target_f) % Default weighting function
			w = get_weights@bt.model.template(self,target_f);
			w(target_f > 25) = 0;
			w(target_f > 8 & target_f < 12) = 100*w(target_f > 8 & target_f < 12);
			w(target_f < 8.5 | target_f > 11) = 0;
			w(target_f > 8.5 & target_f < 11) = 1;
		end

		function [chisq,P] = objective(self,pars)
			P = [];
			chisq = NaN;
			if pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || pars(7)/pars(6) > 20 || any(pars > self.limits(2,:) | pars < self.limits(1,:) )
				return
			else 
				self.set_params(pars);

				% Check stability using is_stable method
				% Will intelligently choose XYZ or Gabcd or Gab stability
				if ~self.p.is_stable()
					return
				end

				% Get nonuniform spectrum
				[~,~,~,~,~,~,~,P] = self.p.spatial_spectrum(self.n_modes,self.n_truncate,self.output_x,self.output_y,self.target_f);

				P = shiftdim(P,2);
				P = bsxfun(@plus,P,1e-9*self.p.emg_a*self.emg*2*pi); % Add the EMG component
				P = bsxfun(@rdivide,P,mex_trapz(self.target_f(self.weights>0),P(self.weights>0,:)));
				P = bsxfun(@times,self.normalization_target,P);
				
				%sqdiff = (abs(P-self.target_P)./self.target_P).^2; % This is the squared fractional difference
				sqdiff = abs(P-self.target_P).^2; % This is the squared fractional difference
				sqdiff = bsxfun(@times,sqdiff,self.weights)/10000;
				chisq = sum(sqdiff(:))/size(P,2); % Failing to divide by the number of electrodes leads to the probability becoming 0 due to numerical precision
			
			end
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
			self.p.apply_variation('t0','cosine_phased',pars(10),pars(11))
		end

		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,mean(target_P,2)); % Fit the first (or only) spectrum given
			a =  db_data.gab(idx,:);
			p = model.params(db_data.iswake(idx));
			initial_values =  [a(1) a(2) a(3)*a(4) a(3)*a(5)*a(7) a(5)*a(8) p.alpha(1) p.beta(1) p.t0 0 0 0.25];
			prior_pp = self.uniform_priors();
		end
		

		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
		end
		
	end
end



