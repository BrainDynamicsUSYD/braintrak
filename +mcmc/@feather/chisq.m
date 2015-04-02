function x = chisq(self,idx)
	self.compress();
	x = [self.fit_data.fitted_chisq];
	if nargin >= 2 && ~isempty(idx)
		x = x(idx);
	end