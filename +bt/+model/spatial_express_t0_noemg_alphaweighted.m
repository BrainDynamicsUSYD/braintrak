classdef spatial_express_t0_noemg_alphaweighted < bt.model.spatial_express_t0_noemg
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
	end

	methods
		function self = spatial_express_t0_noemg_alphaweighted(n_modes) % Constructor
			self.name = 'spatial_express_t0_noemg_alphaweighted';
			self.n_params = 9;
			self.param_names = {'Gee','Gei','Gese','Gesre','Gsrs','Alpha','Beta','t0','t0_amp'};
			self.param_symbols = {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta','t_0','t_0amp'};
			self.param_units = {'','','','','','s^{-1}','s^{-1}','ms','ms'};

			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);
			self.initial_step_size = 0.05.*[2 2 2 2 2 5 10 0.02 0.2];
			self.limits = [
				eps       	20     	  %Gee
				-20     	-eps      %Gei
				eps       	40     	  %Gese
				-40     	-eps      %Gesre
				-14     	-eps      %Gsrs
				10      	100    	  %Alpha
				100     	800    	  %Beta
				75e-3  		140e-3 	  %t0
				-20e-3      20e-3      	  %t0_amp
				].';

			self.p = model.params();
			self.p.disable_set = true;
		end

		function w = get_weights(self,target_f) % Default weighting function
			w = get_weights@bt.model.template(self,target_f);
			w(target_f > 25) = 0;
			w(target_f > 8 & target_f < 12) = 30*w(target_f > 8 & target_f < 12);
		end

		
	end
end



