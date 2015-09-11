classdef reduced_chisqweight < mcmc.model.reduced
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit
	% Note that this does sum over k

	properties
	end


	methods
		function self = reduced_chisqweight() % Constructor
			self = self@mcmc.model.reduced
			self.name = 'model_reduced_chisqweight';
		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			[chisq,P] = objective@mcmc.model.reduced(self,pars);
			if isfinite(chisq)
				sqdiff = (abs(P-self.target_P)./self.target_P).^2;
				chisq = sum(self.target_P.*sqdiff(:).*self.weights(:));
			end
		end
		
	end
end