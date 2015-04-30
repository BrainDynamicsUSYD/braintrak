function f = fit_single(model,f,P,npts,prior_pp,initial_values)
	% This function just uses the default priors to fit a spectrum
	% But it also plots it

	[init_default,prior_default] = model.initialize_fit(f,P);

	if nargin < 6 || isempty(initial_values)
		initial_values = init_default;
	end

	if nargin < 5 || isempty(prior_pp)
		prior_pp = prior_default;
	end

	if nargin < 4 || isempty(npts)
		npts = 24e4;
	end
	
	tic;

	model.set_electrodes('Cz');

	[~,~,f] =  mcmc.fit(model,f,P,prior_pp,initial_values,npts);
	f = mcmc.feather(model,fit_data,plot_data);
	
	fprintf('Fit took %.2f seconds\n',toc);

	f.plot()