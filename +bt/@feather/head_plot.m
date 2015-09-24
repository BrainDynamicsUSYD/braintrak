function head_plot(self,fitted_spatial)
    if ~isa(self.model,'bt.model.template_spatial')
        error('This plot can only be generated for models that inherit from template_spatial')
    end
    
	if nargin < 2
		fitted_spatial = 0;
	end

	if fitted_spatial
		bt_utils.head_imager(self.model.output_x,self.model.output_y,self.fit_data.target_f,self.fit_data.fitted_P,self.model.electrodes,self.fit_data.target_f,self.fit_data.target_P);
	else
		bt_utils.head_imager(self.model.output_x,self.model.output_y,self.fit_data.target_f,self.fit_data.target_P,self.model.electrodes,self.fit_data.target_f,self.fit_data.fitted_P);
	end

