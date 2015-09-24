classdef spatial_emg_xyz < bt.model.template_spatial
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
	end

	methods
		function self = spatial_emg_xyz(n_modes) % Constructor
			self.name = 'spatial_emg_xyz';
			self.n_params = 13;
			self.param_names = {'X','Y','Z','Alpha','Beta','t0','t0_{amp}','EMG_a','EMG_x amp','EMG_x phase','EMG_y amp','EMG_y phase'};
			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.05.*[
				0.2  % X
				0.2  % Y
				0.2  % Z
				40   % Alpha
				100  % Beta
				0.02 % t0
				0.02 % t0_amp
				0.05 % t0_phase
				0.02 % emg_a
				0.01 % emg_x_var
				0.01 % emg_x_phase
				0.01 % emg_y_var
				0.01 % emg_y_phase
			].';
			
			Lx = 0.5;
			head_d = bt.electrode_positions;
			offset = (Lx-head_d)/2;

			self.limits = [
				eps       	1     	 % X
				-1     		1     % Y
				eps       	1     	  % Z
				10      	100    	 % Alpha
				100     	800    	 % Beta
				0.0750  	0.1400 	 % t0
				-0.05           0.05      %t0_amp
				0+offset 	0.5-offset	  % t0_phase
				0 			5		  % emg_a
				-1 			1		  % emg_x_var
				0+offset 	0.5-offset  % emg_x_phase
				-1 			1		  % emg_y_var
				0+offset 	0.5-offset	% emg_y_phase
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
			self.p.apply_variation('t0','cosine_phased',pars(7),pars(8));
			self.p.emg_a = pars(9);
			emg_scaling = 1e-8;
			self.p.apply_variation('emg','cosine_2d_phased',emg_scaling*pars(10),pars(11),emg_scaling*pars(12),pars(13));
		end

		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,squeeze(target_P(:,1))); % Fit the first (or only) spectrum given
			a =  db_data.xyz(idx,:);
			p = model.params(db_data.iswake(idx));
			initial_values =  [a(1) a(2) a(3) p.alpha(1) p.beta(1) p.t0 0 0.25 0 0 0.25 0 0.25];
			prior_pp = self.uniform_priors();
		end	

		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
			xyz(:,1) = params(:,1);
			xyz(:,2) = params(:,2);
			xyz(:,3) = params(:,3);
		end
		
	end
end



