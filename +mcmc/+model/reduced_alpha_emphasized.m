classdef reduced_alpha_emphasized < mcmc.model.reduced
	% This model is the same as reduced, but with increased alpha weighting

	methods
		function self = reduced_alpha_emphasized() % Constructor
			self = self@mcmc.model.reduced
			self.name = 'reduced_alpha_emphasized';
		end

		function w = get_weights(self,target_f) % Default weighting function
			w = get_weights@mcmc.model.reduced(self,target_f);
			w(target_f > 7 & target_f < 12) = 10*w(target_f > 7 & target_f < 12); % Increase alpha weighting
		end
	end
end