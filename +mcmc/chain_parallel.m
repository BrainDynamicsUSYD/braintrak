function [out,chisq_out,accept_ratio] = chain_parallel(model,initial_values,npoints,debugmode,timelimit)
	if nargin < 5 || isempty(timelimit)
		timelimit = false;
	end

	if nargin < 4 || isempty(debugmode)
		debugmode = false;
	end

	% Execute fit_single in parallel according to nworkers and assemble the results into an output matrix

	% Compute the chain
	try
		pool = gcp('nocreate');
		if isempty(pool)
			error('This function should never be be entered into without a pool being open')
		else
			nworkers = pool.NumWorkers;
		end
	catch
		nworkers = matlabpool('size');
		if nworkers == 0
			error('This function should never be be entered into without a pool being open')
		end
	end


	if timelimit
		worker_n = npoints;
	else
		worker_n = round(npoints/nworkers);
		if mod(npoints,nworkers)
			warning('Number of points not divisible by number of workers. Output matrix will be a different size')
		end
	end

	accept = zeros(1,nworkers);

	n_burnin = round(npoints*0.1); % Total burn-in length is based on the final length of the chain
	% This means that if we simulate 1e5 points with however many workers, the BURNIN on each chain 
	% is the same, which means that the result is directly comparable even though the speed-up is sublinear
	% The burnin-length is NOT yet adaptive for parallel fitting yet - proof of concept only

	parfor j = 1:nworkers
		[out1{j},chisq_out1{j},accept(j)] = mcmc.chain(model,initial_values,worker_n,false,timelimit,n_burnin);
	end

	% And assemble
	out = [];
	chisq_out = [];
	for j = 1:nworkers
		out = [out; out1{j}];
		chisq_out = [chisq_out(:); chisq_out1{j}];
	end
	accept_ratio = mean(accept);


	if debugmode
		for j = 1:4
			mcmc.chain_diagnostics(model,out1{j},chisq_out1{j},1,true);
		end
		fprintf(2,'Displaying final chain - note that parallel mode does not show burn-in')
		mcmc.chain_diagnostics(model,out,chisq_out,1,false);
	end
