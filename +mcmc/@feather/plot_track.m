function h_out = plot_track(self,idx,smoothed)
	if nargin < 3 || isempty(idx)
		smoothed = false;
	end

	if nargin < 2 || isempty(idx)
		f = self;
	else
		f = self.subrange(idx);
	end

	if smoothed
		xyz = f.xyz;
		for j = 1:size(xyz,2)
			xyz(:,j) = smooth(xyz(:,j),60);
		end
	else
		xyz = f.xyz;
	end

	h_out = f.plot_statecolored(xyz(:,1),xyz(:,2),xyz(:,3));
