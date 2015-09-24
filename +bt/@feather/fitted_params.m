function x = fitted_params(self)
	self.compress();
	x = reshape([self.fit_data.fitted_params],self.model.n_params,[]).';