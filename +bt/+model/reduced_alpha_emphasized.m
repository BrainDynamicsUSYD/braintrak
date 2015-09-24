classdef reduced_alpha_emphasized < bt.model.reduced
	% This model is the same as reduced, but with increased alpha weighting

	methods
		function self = reduced_alpha_emphasized() % Constructor
			self = self@bt.model.reduced
			self.name = 'reduced_alpha_emphasized';
		end

		function w = get_weights(self,target_f) % Default weighting function
			w = get_weights@bt.model.reduced(self,target_f);
			w(target_f > 6 & target_f < 10) = 25*w(target_f > 6 & target_f < 10); % Increase alpha weighting
		end
	end
end