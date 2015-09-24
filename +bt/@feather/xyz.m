function x = xyz(self)
	self.compress();
	x = reshape([self.fit_data.xyz],3,[]).';