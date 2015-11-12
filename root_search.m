function [f_out,fit_data,plot_data] = root_search(model,npoints,debugmode)
	% Take in a model, spectrum, prior distribution, initial values, and chain length
	% Return a feather object
	% The posterior distribution correctly reflects the skip_fit state of the model
	% This function is intended to be called via a wrapper like fit_single
		
	% Correctly set the skip_fit state of the model and normalization target
	% As a hacky workaround - can supply npoints as '15s' for a time limit
	tic;

	if nargin < 3 || isempty(debugmode)
		debugmode = false;
	end

	if nargin < 2 || isempty(npoints)
		npoints = '10s'; % Default run time
	end

	if ischar(npoints)
		timelimit = sscanf(npoints,'%fs');
		npoints = 1e5; % Default burning - this needs to be tuned
	else
		timelimit = false;
	end

	target_f = 1;
	target_P = 1;
	skip_fit = [];

	[initial_values,prior_pp] = model.initialize_fit(target_f,target_P);

	if isempty(model.electrodes)
		model.set_electrodes('Cz');
	end

	model.prepare_for_fit(target_f,target_P,initial_values,prior_pp,skip_fit);

	try
		pool = gcp('nocreate');
	catch
		nworkers = matlabpool('size');
		if nworkers == 0
			pool = [];
		else
			pool.NumWorkers = nworkers;
		end
	end

	if isempty(pool)
		fprintf('Single CPU mode\n');
		[out,posterior_out,accept_ratio] = bt.core.chain(model,initial_values,npoints,debugmode,timelimit);
	else
		fprintf('Parallel mode with %d workers\n',pool.NumWorkers);
		[out,posterior_out,accept_ratio] = bt.core.chain_parallel(model,initial_values,npoints,debugmode,timelimit);
	end

	nroots = 5.1+log(posterior_out)/2;

	fprintf('Tested %d points\n',length(posterior_out))
	a =  histc(nroots,0:5);
	for j = 0:5
		fprintf('%d root(s) : %d\n',j,a(j+1));
	end

	keyboard

