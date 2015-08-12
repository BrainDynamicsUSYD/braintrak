function h_out = plot_track(self,idx,smoothed)
	if nargin < 3 || isempty(idx)
		smoothed = 1; % With a moving average of 1, smoothing does nothing
	end

	if nargin < 2 || isempty(idx)
		idx = 1:self.latest;
	end

	% First, do the smoothing
	f = self.subrange(idx);
	xyz = f.xyz;
	for j = 1:size(xyz,2)
		xyz(:,j) = smooth(xyz(:,j),smoothed);
	end

	% Now skip the affected parts
	output_indexes = idx(ceil(smoothed/2):(end-(ceil(smoothed/2)-1)));
	fprintf('Skipping %d samples - total length = %d samples',length(idx)-length(output_indexes),length(output_indexes))
	if isempty(output_indexes)
		error('Too much smoothing for provided indexes');
	end



	f2 = self.subrange(output_indexes);
	xyz = xyz(ceil(smoothed/2):(end-(ceil(smoothed/2)-1)),:);
	h_out = f2.plot_statecolored(xyz(:,1),xyz(:,2),xyz(:,3));
