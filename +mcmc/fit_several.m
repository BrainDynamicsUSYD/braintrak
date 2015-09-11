function fdata = fit_several(model,dataset,subject_idx,t_increment,npts_per_fit,use_prior)
	% CALLING SEQUENCE
	% fit_subject loads subject data
	% fit_spectrum takes in a spectrum and returns a fitted point struct
	% fit_single(_parallel) actually does the fitting of the 8 parameters

	if nargin < 6 || isempty(use_prior)
		use_prior = true;
	end

	debugmode = false;

	if isa(model,'mcmc.model.spatial_t0_2d') || isa(model,'mcmc.model.template_spatial')
		d = mcmc.load_subject_data(dataset,subject_idx,{'all'});
		disp('Multi-electrode fitting selected')
	else
		d = mcmc.load_subject_data(dataset,subject_idx,{'Cz'});
		disp('Single-electrode fitting selected')
	end

	model.set_electrodes(d);

	if nargin < 5 || isempty(npts_per_fit)
		npts_per_fit = 1e5;
		fprintf(1,'Default number of points: %d\n',npts_per_fit)
	end

	if nargin < 4 || isempty(t_increment)
		t_increment = 1:1:5; % 3600=1 hour. This is an array of offsets
	end
	
	initial_idx = 1;

	for j = 1:length(t_increment)
		tic;
		fprintf('Fitting: %d/%d\n',j,length(t_increment))
		target_state{j} = d.state_str{initial_idx+(t_increment(j)-1)};
		target_P = squeeze(d.s(:,initial_idx+(t_increment(j)-1),:));

		if j == 1
			[initial_params,initial_pp] = model.initialize_fit(d.f,target_P);
			[fit_data(j),plot_data(j)] = mcmc.fit(model,d.f(:),target_P,initial_pp,initial_params,npts_per_fit,target_state{j},[],debugmode);
			fdata = mcmc.feather(model,fit_data(j),plot_data(j),t_increment(j));
		else
			if use_prior
				[fit_data(j),plot_data(j)] = mcmc.fit(model,d.f(:),target_P,fit_data(j-1).posterior_pp,fit_data(j-1).fitted_params,npts_per_fit,target_state{j},[],debugmode);			
			else
				[fit_data(j),plot_data(j)] = mcmc.fit(model,d.f(:),target_P,initial_pp,fit_data(j-1).fitted_params,npts_per_fit,target_state{j},[],debugmode);			
            end
            fdata.insert(fit_data(j),plot_data(j),t_increment(j));
		end
	%	fdata.plot(j)
	end
	toc;
	fdata.compress()