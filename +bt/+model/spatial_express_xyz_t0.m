classdef spatial_express_xyz_t0 < bt.model.template_spatial
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
	end

	methods
		function self = spatial_express_xyz_t0(n_modes) % Constructor
			self.name = 'spatial_express_xyz_t0';
			self.n_params = 8;
			self.param_names = {'X','Y','Z','Alpha','Beta','t0','EMGa','t0_{amp}'};
			self.param_symbols = {'X','Y','Z','\alpha','\beta','t_0','A_{EMG}','t_{0amp}'};
			self.param_units = {'','','','s^{-1}','s^{-1}','ms','','ms'};

			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.05.*[0.2 0.2 0.2 5 10 0.02 0.05 0.02];
			self.limits = [
				eps       	1     	 
				-1     		1     
				eps       	1     	  
				10      	100    	 
				100     	800    	 
				0.0750  	0.1400 	 
				0       	5      	 
				-0.05       0.05      	  %t0_amp
				].';

			self.p = model.params();
			self.p.disable_set = true;
		end

		function set_params(self,pars)
			% This function is what decides what the spatial variations are actually going to be
			%e.g. p.apply_variation('cosine','t0',pars(10));
			self.p.xyz = pars(1:3);
			self.p.alpha(:) = pars(4);
			self.p.beta(:) = pars(5);
			self.p.t0 = pars(6);
			self.p.taues = self.p.t0/2;
			self.p.tause = self.p.t0/2;
			self.p.emg_a = pars(7);
			self.p.apply_variation('t0','cosine',pars(8));
		end

		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,squeeze(target_P(:,1))); % Fit the first (or only) spectrum given
			a =  db_data.xyz(idx,:);
			p = model.params(db_data.iswake(idx));
			initial_values =  [a(1) a(2) a(3) p.alpha(1) p.beta(1) p.t0 0 0];
			prior_pp = self.uniform_priors();
		end	

		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1);
			xyz(:,2) = params(:,2);
			xyz(:,3) = params(:,3);
		end
		
	end
end



