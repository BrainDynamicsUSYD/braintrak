function x = state_str(self)
	self.compress();
	x = {self.fit_data.state_str};