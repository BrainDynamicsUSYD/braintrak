function [f,bic] = fit_br_bic(idx)
	romesh_utils.matlabpool_cluster(12)
	tic;
	eegdb = data.eeg_database;
	files = eegdb.get('br_tfs','EC');
	eegdb.close();

	n_fits = length(files);
	bic = zeros(n_fits,1);
	f(n_fits) = bt.feather;

	debugmode = false;
	npts_per_fit = 5e4;

	parfor j = 1:n_fits

		initial_params_override = [];
		m = [];
		skip_fit = [];

		switch idx
			case 1 % Full model
				m = bt.model.full;
				skip_fit = [0 0 0 0 0 0 0 0 0];
			case 2
				m = bt.model.full;
				skip_fit = [0 0 0 0 0 0 0 0 2];
				initial_params_override = [NaN NaN NaN NaN NaN NaN NaN NaN eps ];
			case 3
				m = bt.model.full_emgf;
				skip_fit = [0 0 0 0 0 0 0 0 0 0];
			case 4 
				m = bt.model.reduced;
				skip_fit = [0 0 0 0 0 0 0];
			case 5
				m = bt.model.reduced; 
				skip_fit = [0 0 0 0 0 0 2];
				initial_params_override = [NaN NaN NaN NaN NaN NaN eps ];
			case 6
				m = bt.model.reduced_emgf; 
				skip_fit = [0 0 0 0 0 0 0 0];
			otherwise
				error('Undefined index')
		end

		d = load(files{j});
		m.set_electrodes(d.electrode_choice);
		[initial_params,initial_pp] = m.initialize_fit(d.f(:),d.P(:));
		if ~isempty(initial_params_override)
			initial_params(isfinite(initial_params_override)) = initial_params_override(isfinite(initial_params_override));
		end

		f(j) = bt.fit_spectrum(m,d.f(:),d.P(:),initial_pp,initial_params,npts_per_fit,[],skip_fit,debugmode);

		f(j).compress();
		bic(j) = f(j).fit_data.bic;
	end

	save(sprintf('bic_master_case_%d',idx),'bic','f','-v7.3');
	toc;

