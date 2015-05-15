function h_out = plot_track(self,idx,smoothed)
	if nargin < 3 || isempty(idx)
		smoothed = false;
	end

	if nargin < 2 || isempty(idx)
		f = self;
	else
		f = self.subrange(idx);
	end

	smoothed = true;
	if smoothed
		% pars = f.fitted_params;
		% for j = 1:size(pars,2)
		% 	pars(:,j) = smooth(pars(:,j),60);
		% end
		% xyz = f.model.get_xyz(pars);

		xyz = f.xyz;
		for j = 1:size(xyz,2)
			xyz(:,j) = smooth(xyz(:,j),60);
		end
	else
		xyz = f.xyz;
	end

	contaminated = braintrack_utils.chisq_outliers(f.chisq);
	xyz(contaminated,:) = NaN;
	h_out = f.plot_statecolored(xyz(:,1),xyz(:,2),xyz(:,3));
