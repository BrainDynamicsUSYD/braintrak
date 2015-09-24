function params = fitted_params_from_posterior(self,xyz_mode)
	% feather.fitted_params returns the maximum likelihood fitted parameters
	% This function returns the parameter value corresponding to the maximum 
	% value of the posterior density. This value reflects the entire distribution
	% and should be more stable (i.e. smoother). However, if these parameters are
	% used directly in spectrum, the result will probably NOT be as good as the
	% actual fitted parameters e.g. due to interdependencies. 
	% For XYZ mode, use fit_data.posterior_xyz instead of fit_data.posterior_pp

	if nargin < 2 || isempty(xyz_mode)
		xyz_mode = false;
	end

	self.compress();

	if xyz_mode
		n = 3;
	else
		n = self.model.n_params;
	end

	params = zeros(self.latest,n);

	for j = 1:self.latest
		if xyz_mode
			post = self.fit_data(j).xyz_posterior;
		else
			post = self.fit_data(j).posterior_pp;
		end

		%for k = 1:self.model.n_params
			% Could fit the curve here instead
			%gauss = @(pars,x) pars(1).*exp((-(x-pars(2)).^2)/(2*pars(3).^2));
		%end
		[~,idx] = max(post.y);
		params(j,:) =  arrayfun(@(j) post.x(idx(j),j),1:n);
	end