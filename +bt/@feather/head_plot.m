function head_plot(self,fitted_spatial)
	if nargin < 2
		fitted_spatial = 0;
	end

	if fitted_spatial
		multielectrode.head_imager(self.model.output_x,self.model.output_y,self.fit_data.target_f,self.fit_data.fitted_P,self.model.electrodes,self.fit_data.target_f,self.fit_data.target_P);
	else
		multielectrode.head_imager(self.model.output_x,self.model.output_y,self.fit_data.target_f,self.fit_data.target_P,self.model.electrodes,self.fit_data.target_f,self.fit_data.fitted_P);
	end

