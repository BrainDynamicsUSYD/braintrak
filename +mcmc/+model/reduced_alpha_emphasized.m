classdef reduced_alpha_emphasized < mcmc.model.reduced
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit
	% Note that this does sum over k

	properties
	end


	methods
		function self = reduced_alpha_emphasized() % Constructor
			self = self@mcmc.model.reduced
			self.name = 'reduced_alpha_emphasized';
		end

		function w = get_weights(self,target_f) % Default weighting function
			w = ones(size(target_f));
			w = target_f.^-1; % Decades equally weighted
			w(target_f > 8 & target_f < 12) = 15*w(target_f > 8 & target_f < 12); % Increased alpha weighting
		end
	end
end