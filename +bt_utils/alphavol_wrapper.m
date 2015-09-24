function h = alphavol_wrapper(a,r)
	if nargin < 2 || isempty(r)
		r = Inf;
	end

	[V,S] = utils.alphavol(a,r);
	k = S.bnd;
	map = unique(k(:));
	map(:,2) = 1:length(map); % First col is the index in k, second row is the new index
	vertices = a(map(:,1),:); % A final list of vertices
	faces = arrayfun(@(b) map(find(map(:,1)==b),2),k);
	h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
