function fdata = fit_kaggle(model,set_type,set_number,npts_per_fit)
	% CALLING SEQUENCE
	% fit_subject loads subject data
	% fit_spectrum takes in a spectrum and returns a fitted point struct
	% fit_single(_parallel) actually does the fitting of the 8 parameters

	%system(sprintf('mkdir -p %s_testfit2',subject));
	debugmode = false;
	matlabpool_cluster(12)
	%npts_per_fit = '10s';

	% if isa(model,'bt.model.spatial_t0_2d') || isa(model,'bt.model.template_spatial')
	% 	d = bt.core.load_subject(dataset,subject_idx,{'all'});
	% 	disp('Multi-electrode fitting selected')
	% else
	% 	d = bt.core.load_subject(dataset,subject_idx,{'Cz'});
	% 	disp('Single-electrode fitting selected')
	% end
	d = bt.load_kaggle_data(set_type,set_number);

	model.set_electrodes(d);
	model.k0 = Inf; % Disable volume conduction for intracranial recordings

	t_increment = 1:1:3569; % Fit the recording minutes of the recording
	initial_idx = 1;

	for j = 1:length(t_increment)
		tic;
		fprintf('Fitting: %d/%d\n',j,length(t_increment))
		target_state{j} = d.state_str{initial_idx+(t_increment(j)-1)};
		target_P = squeeze(d.s(:,initial_idx+(t_increment(j)-1),:));

		if j == 1
			[initial_params,initial_pp] = model.initialize_fit(d.f,target_P);
			[~,fit_data(j),plot_data(j)] = bt.core.fit_spectrum(model,d.f(:),target_P,initial_pp,initial_params,npts_per_fit,target_state{j},[],debugmode);
			fdata = bt.feather(model,fit_data(j),plot_data(j),t_increment(j));
		else
			[~,fit_data(j),plot_data(j)] = bt.core.fit_spectrum(model,d.f(:),target_P,fit_data(j-1).posterior_pp,fit_data(j-1).fitted_params,npts_per_fit,target_state{j},[],debugmode);			
			%[fit_data(j),plot_data(j)] = bt.core.fit_spectrum(model,d.f(:),target_P,initial_pp,fit_data(j-1).fitted_params,npts_per_fit,target_state{j},[],debugmode);			
            fdata.insert(fit_data(j),plot_data(j),t_increment(j));
		end

		if ~mod(j,30)
			save(sprintf('kaggle_temp_%s_%d',set_type,set_number),'fdata');
		end
		%fdata.plot(j)
		%drawnow
	end
	toc;
	fdata.compress();
	f = fdata;
	save(sprintf('kaggle_%s_%d',set_type,set_number),'f');
