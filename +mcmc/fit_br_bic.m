function [f,bic] = fit_br_bic(idx)
	romesh_utils.matlabpool_cluster(12)
	tic;
	eegdb = data.eeg_database;
	files = eegdb.get('br_tfs','EC');
	eegdb.close();

	n_fits = length(files);
	bic = zeros(n_fits,1);
	f(n_fits) = mcmc.feather;

	debugmode = false;
	npts_per_fit = 1e5;

	parfor j = 1:n_fits

		initial_params_override = [];
		m = [];
		skip_fit = [];

		switch idx
			case 1 % Full model
				m = mcmc.model.full;
				skip_fit = [0 0 0 0 0 0 0 0 0];
			case 2
				m = mcmc.model.full;
				skip_fit = [0 0 0 0 0 0 0 0 2];
				initial_params_override = [NaN NaN NaN NaN NaN NaN NaN NaN 0 ];
			case 3
				m = mcmc.model.full_b_ratio;
				skip_fit = [0 0 0 0 0 0 2 0 0];
				initial_params_override = [NaN NaN NaN NaN NaN NaN 4 NaN NaN ];
			case 4
				m = mcmc.model.full_b_ratio;
				skip_fit = [0 0 0 0 0 0 2 0 0];
				initial_params_override = [NaN NaN NaN NaN NaN NaN 10 NaN NaN ];
			case 5
				m = mcmc.model.full_b_ratio;
				skip_fit = [0 0 0 0 0 0 2 0 2];
				initial_params_override = [NaN NaN NaN NaN NaN NaN 4 NaN 0 ];
			case 6
				m = mcmc.model.full_b_ratio;
				skip_fit = [0 0 0 0 0 0 2 0 2];
				initial_params_override = [NaN NaN NaN NaN NaN NaN 10 NaN 0 ];
			case 7 
				m = mcmc.model.reduced;
				skip_fit = [0 0 0 0 0 0 0];
			case 8
				m = mcmc.model.reduced; 
				skip_fit = [0 0 0 0 0 0 2];
				initial_params_override = [NaN NaN NaN NaN NaN NaN 0 ];
			case 9
				m = mcmc.model.reduced_b_ratio; 
				skip_fit = [0 0 0 0 2 0 0];
				initial_params_override = [NaN NaN NaN NaN 4 NaN NaN ];
			case 10
				m = mcmc.model.reduced_b_ratio; 
				skip_fit = [0 0 0 0 2 0 0];
				initial_params_override = [NaN NaN NaN NaN 10 NaN NaN ];			
			case 11 
				m = mcmc.model.reduced_b_ratio; 
				skip_fit = [0 0 0 0 2 0 2];
				initial_params_override = [NaN NaN NaN NaN 4 NaN 0 ];
			case 12 
				m = mcmc.model.reduced_b_ratio; 
				skip_fit = [0 0 0 0 2 0 2];
				initial_params_override = [NaN NaN NaN NaN 10 NaN 0 ];
			otherwise
				error('Undefined index')
		end

		d = load(files{j});
		m.set_electrodes(d.electrode_choice);
		[initial_params,initial_pp] = m.initialize_fit(d.f(:),d.P(:));
		if ~isempty(initial_params_override)
			initial_params(isfinite(initial_params_override)) = initial_params_override(isfinite(initial_params_override));
		end
		
		[~,~,f(j)] = mcmc.fit(m,d.f(:),d.P(:),initial_pp,initial_params,npts_per_fit,[],skip_fit,debugmode);
		f(j).compress();
		bic(j) = f(j).fit_data.bic;
	end

	save(sprintf('bic_lf_case_%d_all',idx),'bic','f','-v7.3');
	toc;

