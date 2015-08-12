function f = fit_br_spatial(model,subject_id,npts)
	fdata = mcmc.load_br_data(subject_id,'all');

	[init_default,prior_pp] = model.initialize_fit(fdata.f,fdata.P);

	model.set_electrodes(fdata.electrodes);
	[~,~,f] =  mcmc.fit(model,fdata.f,fdata.P,prior_pp,init_default,npts,[],[]);

	fprintf('Fit took %.2f seconds\n',toc);

	%f.plot()