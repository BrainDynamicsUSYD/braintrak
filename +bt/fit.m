function [fit_data,plot_data,f_out] = fit(model,target_f,target_P,prior_pp,initial_values,npoints,target_state,skip_fit,debugmode)
	% Take in a model, spectrum, prior distribution, initial values, and chain length
	% Return the fitted parameters and posterior distribution
	% The posterior distribution correctly reflects the skip_fit state of the model
	% This function is intended to be called via a wrapper like fit_single
		
	% Correctly set the skip_fit state of the model and normalization target
	% As a hacky workaround - can supply npoints as '15s' for a time limit
	tic;

	if nargin < 9 || isempty(debugmode)
		debugmode = false;
	end

	if nargin < 8 || isempty(skip_fit)
		skip_fit = [];
	end

	if nargin < 7
		target_state = 'N/A';
	end

	if ischar(npoints)
		timelimit = sscanf(npoints,'%fs');
		npoints = 1e5; % Default burning - this needs to be tuned
	else
		timelimit = false;
	end

	assert(size(target_P,2)==length(model.electrodes))
	
	model.prepare_for_fit(target_f,target_P,initial_values,prior_pp,skip_fit);

	% Compute the chain
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
		[out,posterior_out,accept_ratio] = bt.chain(model,initial_values,npoints,debugmode,timelimit);
	else
		fprintf('Parallel mode with %d workers\n',pool.NumWorkers);
		[out,posterior_out,accept_ratio] = bt.chain_parallel(model,initial_values,npoints,debugmode,timelimit);
	end

	% Pick and evaluate the fitted parameters
	fit_data.target_f = target_f;
	fit_data.target_P = target_P;
	fit_data.skip_fit = model.skip_fit;
	fit_data.chain_length = length(posterior_out);
	[~,b] = max(posterior_out);
	fit_data.fitted_params = out(b,:);
	fit_data.fitted_posterior = posterior_out(b);
	fit_data.posterior_pp = model.make_posterior(out);
	fit_data.xyz_posterior = model.xyz_posterior(out);
	
	[fit_data.fitted_chisq,~,likelihood,fit_data.fitted_P] = model.probability(fit_data.fitted_params);
	
	fit_data.bic = -2*log(likelihood)+model.n_fitted*log(sum(model.weights>0));
	fit_data.aic = 2*model.n_fitted - 2*log(likelihood);
	fit_data.aicc = fit_data.aic + (2*model.n_fitted*(model.n_fitted+1))/(sum(model.weights>0) - model.n_fitted - 1);

	fit_data.xyz = model.get_xyz(fit_data.fitted_params);
	fit_data.fit_time = toc;
	fit_data.state_str = target_state;

	% If the model is calling for the posterior distribution to be held constant
	% then copy the prior distribution into the posterior distribution 
	for j = 1:length(fit_data.skip_fit)
		if fit_data.skip_fit(j)
			fit_data.posterior_pp.x(:,j) = prior_pp.x(:,j);
			fit_data.posterior_pp.y(:,j) = prior_pp.y(:,j);
			fit_data.posterior_pp.ndx(:,j) = prior_pp.ndx(:,j);
		end
	end

	% Compute the Vcount matrix for plotting later
	plot_data.xyzlim = [0  -1     0 ; 1   1   1.3];
	plot_data.gridres = 0.025;
	xyz = model.get_xyz(out);
	plot_data.Vcount = utils.simple_bin3d(xyz,plot_data.xyzlim,plot_data.gridres);
	plot_data.accept_ratio = accept_ratio;

	p = model.p_from_params(fit_data.fitted_params);
	[plot_data.tent_x,plot_data.tent_y,plot_data.tent_z,plot_data.tent_u] = tent.compute(p.alpha(1),p.beta(1),p.t0,p.gammae);
	plot_data.tent_alpha = +(plot_data.tent_x+plot_data.tent_y<0.90 & plot_data.tent_x > 0 & plot_data.tent_x < 1);
	
	if nargout >= 3
		f_out = bt.feather(model,fit_data,plot_data);
	end

