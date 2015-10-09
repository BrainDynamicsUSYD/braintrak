function [keep,output_params] = nd_cutdown(params,gridres,npoints,lowlim)
	% Perform uniform resampling in bins
	% Params- matrix of parameter values
	% gridres- bin size in parameter units
	% npoints- number of points per bin

	%params = randn(100000,8);

	%gridres = 0.20;
	%npoints = 5;

	% First, need to calculate the LOWER limit
	if nargin < 4 || isempty(lowlim)
		lowlim = floor(min(params)/gridres)*gridres; % This is the lower limit on the bins
	end

	voxels = zeros(size(params));
	for j = 1:size(params,2) % For each parameter, calculate the bin
		voxels(:,j) = floor((params(:,j)-lowlim(j))/gridres)+1;
	end

	[voxels,sort_order] = sortrows(voxels);
	[~,~,check] = unique(voxels,'rows'); % These are only the voxels which need checking

	if nargin < 2 || isempty(npoints)
		% automatically set npoints so that there are about 1% of the original points
		npoints = floor(size(params,1)/100/length(unique(check)));
	end

	rng = RandStream('mt19937ar','Seed',10);
	keep = logical(zeros(1,size(params,1)));

	%check(end+1) = max(check)+1; % Add a dummy entry at the end
	stop_idx = find(diff([check ; max(check+1)])); % These are the stop indices
	start_idx = [1 ; stop_idx(1:end-1)+1];

	for j = 1:length(start_idx) % for each voxel
		matching = start_idx(j):stop_idx(j);
		
		if length(matching) > npoints % Randomly select points from bin
			matching = randsample(rng,matching,npoints);
		end
		keep(matching) = 1;
	end
	[~,lookup] = sort(sort_order);
	keep = keep(lookup);

	fprintf(1,'Selected %d points of %d (%.1f%%)\n',sum(keep),length(keep),sum(keep)/length(keep)*100);
	
	if nargout > 1
		output_params = params(keep,:);
	end

	%return

	% figure
	% analysis.cloud_plots(params(:,1),params(:,2),[],[],[100 100])
	
	% figure
	% analysis.cloud_plots(params(keep,1),params(keep,2),[],[],[100 100])

	% figure
	% [idx,dist] = knnsearch(params,params,'k',2);
	% ksdensity(dist(:,2))
	% title('original')
	
	% figure
	% [idx,dist] = knnsearch(params(keep,:),params(keep,:),'k',2);
	% ksdensity(dist(:,2))
	% title('cutdown')

	% keyboard