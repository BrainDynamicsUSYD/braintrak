classdef reduced_ly_no_n < bt.model.reduced_ly
	% This model is set up to fit t0 keeping the prior fixed if the alpha peak is not present
	% Note that prepare_for_fit must be used before spectrum() can be used
	% Uniform priors and DB initial fit
	% Note that this does sum over k

	properties
	end


	methods
		function self = reduced_ly_no_n() % Constructor
			self = self@bt.model.reduced_ly
			self.name = 'model_reduced_ly_no_n';
		end

		function [chisq,P] = objective(self,pars) % Calculate the objective
			chisq = NaN;
			P = [];
			if any(pars > self.limits(2,:) | pars < self.limits(1,:))
				return
			end

		    if self.compute_Mtot % If this flag has been set to false, then these fields are already present
		        % Should be fine to let that just throw an undefined variable error later if this isn't the case
		       self.stab_Mtot = exp(1i*self.stab_w*pars(6));
		       self.spec_Mtot = exp(1i*self.spec_w*pars(6));
		    end

			stab_L = 1./((1-1i*self.stab_w/pars(4)).*(1-1i*self.stab_w/pars(5)));
		    spec_L = 1./((1-1i*self.spec_w/pars(4)).*(1-1i*self.spec_w/pars(5)));
		    zprime = (pars(4)+pars(5))^2/(pars(4)*pars(5))*pars(3); % Effectively Gsrs

		    % cf. James corticothalamic dynamics eq 9
		    d = ((self.stab_gamma_prefactor) - pars(1)).*(1+stab_L.*stab_L.*zprime) - stab_L.*stab_L.*pars(2).*(1+zprime).*self.stab_Mtot;
		    stab=d(1)>0 && ~any(real(d(2:end))<0 & imag(d(2:end)).*imag(d(1:end-1))<0);

		    if ~stab
		    	return
		    end

		    % And the spectrum. At k = 0
		    one_plus_zprime = (1+zprime.*spec_L.*spec_L);

		    P = zeros(size(self.spec_w));
		    q2re2 = (self.spec_gamma_prefactor - pars(1) -(spec_L.*spec_L.*pars(2)*(1+zprime))./one_plus_zprime.*self.spec_Mtot);
		    k2u = self.k2u;
		    k2_volconduct = self.k2_volconduct;
		    re2 = self.p.re.^2;

		    for j = 1:size(k2u,1)
		        P = P + k2u(j,2).*abs((1)./(one_plus_zprime.*(k2u(j,1)*re2+q2re2))).^2 * k2_volconduct(j); % For use with the efficient way
		    end

		    P = P + 0*pars(7)*self.emg; % Add the EMG component
		    
		    P = P./mex_trapz(self.target_f(self.weights>0),P(self.weights>0));
		    P = self.normalization_target*P;
		    sqdiff = (abs(P-self.target_P)./self.target_P).^2; % This is the squared fractional difference
		    chisq = sum(sqdiff(:).*self.weights(:));
		end
		
	end
end