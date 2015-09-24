function fdata = speed_test(model,npts_per_fit)
	% CALLING SEQUENCE
	% fit_subject loads subject data
	% fit_spectrum takes in a spectrum and returns a fitted point struct
	% fit_single(_parallel) actually does the fitting of the 8 parameters

	%system(sprintf('mkdir -p %s_testfit2',subject));
	debugmode = false;
	%npts_per_fit = 1e5;

	d = load('seizure_test.mat');
	padding_minutes = 5; % Start this many minutes before sleep onset
	start_idx = max(1,find(d.state_score>2,1,'first')-padding_minutes*60);

	d.t = d.t(start_idx:end);
	d.nspec = d.nspec(start_idx:end,:);
	d.s = d.s(:,start_idx:end,:);

	d.s = bsxfun(@rdivide,d.s,trapz(d.f,d.s));
	

	d.state_str = d.state_str(start_idx:end);
	d.state_score = d.state_score(start_idx:end);
	d.start_idx = start_idx;

	model.set_electrodes('Cz');
	max_iters = 20; % 3600=1 hour
	initial_idx = 1;

	figure
	e_plot = loglog(d.f,d.s(:,1));
	hold on
	f_plot = loglog(d.f,d.s(:,1),'r--');
	set(gca,'XLim',[1 45])
	
	for j = 1:max_iters
		tic;
		fprintf('Fitting: %d/%d\n',j,max_iters)
		
		target_state{j} = d.state_str{initial_idx+(j-1)};
		target_P = squeeze(d.s(:,initial_idx+(j-1)));

		if j == 1
			[initial_params,initial_pp] = model.initialize_fit(d.f,target_P);
			[fit_data(j),plot_data(j)] = bt.fit(model,d.f(:),target_P,initial_pp,initial_params,npts_per_fit,target_state{j},[],debugmode);
			fdata = bt.feather(model,fit_data(j),plot_data(j));
		else
			[fit_data(j),plot_data(j)] = bt.fit(model,d.f(:),target_P,fit_data(j-1).posterior_pp,fit_data(j-1).fitted_params,npts_per_fit,target_state{j},[],debugmode);			
            fdata.insert(fit_data(j),plot_data(j));
		end

		set(e_plot,'YData',target_P);
		set(f_plot,'YData',fdata.fit_data(fdata.latest).fitted_P);
		drawnow
	end
	toc;
	fdata.compress()
	fdata.plot
