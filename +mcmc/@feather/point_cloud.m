function point_cloud(self)
	% Plots a color-coded point cloud
	xyz = self.xyz;
	figure
	self.scatter_statecolored(xyz(:,1),xyz(:,2),xyz(:,3));
	tent.surface_hg2([],false,'a1')