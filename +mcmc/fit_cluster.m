function fdata = fit_cluster(model,dataset,subject_idx,npts_per_fit,suffix)
	% if matlabpool('size') == 0
	% 	matlabpool_cluster(12)
	% end

	debugmode = false;
	%npts_per_fit = 12e5;


	output_dir = sprintf('%s_%i_%s',dataset,subject_idx,suffix);
	if exist(sprintf('%s/fit_output.mat',output_dir)) % Resume an interrupted run
		error('A finalized fit already exists!');
	end

	if isa(model,'mcmc.model.spatial_t0_2d') || isa(model,'mcmc.model.template_spatial')
		d = mcmc.load_subject_data(dataset,subject_idx,{'all'});
		disp('Multi-electrode fitting selected')
	else
		d = mcmc.load_subject_data(dataset,subject_idx);
		disp('Single-electrode fitting selected')
	end

	model.set_electrodes(d);

	t_increment = 1:1:length(d.t); % This is the most we can do
	% Find the first increment with a clean spectrum
	clean_idx = find(d.nspec >= 20,1,'first');
	if isempty(clean_idx) || clean_idx > t_increment(end)
		error('No clean spectrum to start');
	else
		t_increment(t_increment<clean_idx) = [];
	end

	initial_idx = 1;

	% Handle the caching
	if exist(sprintf('%s/cached_output.mat',output_dir)) % Resume an interrupted run
		% All of the arguments will be overwritten by this step, as they should be
		j = 0; % Pre-initialize
 		try
			load(sprintf('%s/cached_output.mat',output_dir));
		catch
			load(sprintf('%s/cached_output_backup.mat',output_dir));
		end
		j_start = j+1;
	else % Initialize variables
		mkdir(output_dir);
		mkdir(sprintf('%s/plot_data',output_dir));

		j_start = 1;
		target_state{j_start} = d.state_str{j_start+(t_increment(j_start)-1)};
		target_P = squeeze(d.s(:,j_start+(t_increment(j_start)-1),:));
		[initial_params,initial_pp] = model.initialize_fit(d.f,target_P);
		[fit_data,plot_data] = mcmc.fit(model,d.f(:),target_P,initial_pp,initial_params,npts_per_fit,target_state{initial_idx},[],debugmode);
		plot_data_fname = sprintf('%s/plot_data/t_%d',output_dir,j_start);
		save(plot_data_fname,'plot_data');
		f = mcmc.feather(model,fit_data,plot_data_fname);
		j_start = 2;
	end

	tic;
	for j = j_start:length(t_increment)
		fprintf('Fitting: %d/%d\n',j,length(t_increment))
		didx = initial_idx+(t_increment(j)-1); % Current index of d
		target_P = squeeze(d.s(:,didx,:));

		if d.nspec(didx) < 20 % If there were more than 10s marked bad in the data
			fit_data.skip_fit(1:end) = 3;
			fit_data.target_f = d.f(:);
			fit_data.target_P = target_P;
			fprintf('Skipping %i due to bad data\n',j);
		else
			[fit_data,plot_data] = mcmc.fit(model,d.f(:),target_P,fit_data.posterior_pp,fit_data.fitted_params,npts_per_fit,d.state_str{didx});	
			fprintf('Fitted %i\n',j);
		end

		plot_data_fname = sprintf('%s/plot_data/t_%d',output_dir,j);
		save(plot_data_fname,'plot_data');
		f.insert(fit_data,plot_data_fname,t_increment(j)); % Insert after saving as j will be recomputed

		if ~mod(j,3) % Save backup after every 30s chunk
			fprintf('Saving backup: Just fitted t=%ds\n',t_increment(j))
			backup_name = sprintf('%s/cached_output.mat',output_dir);
			system(sprintf('cp %s %s',backup_name,strrep(backup_name,'cached_output','cached_output_backup')));
			ts1 = toc;
			f.compress();
			save(sprintf('%s/cached_output',output_dir),'f','j','fit_data','plot_data','t_increment');
			ts2 = toc;
			fprintf('Finished saving: Took %.2fs\n',ts2-ts1);
		end
	end
	toc;
	f.compress()
	save(sprintf('%s/fit_output',output_dir),'f');
