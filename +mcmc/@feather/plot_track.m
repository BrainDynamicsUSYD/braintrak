function h_out = plot_track(self,idx)
	if nargin < 2 || isempty(idx)
		f = self;
	else
		f = self.subrange(idx);
	end

	xyz = f.xyz;
	contaminated = tracking.chisq_outliers(f.chisq);
	xyz(contaminated,:) = NaN;
	h_out = self.plot_statecolored(xyz(:,1),xyz(:,2),xyz(:,3));
