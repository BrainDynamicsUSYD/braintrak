function fdata = fit_seizure(subject_idx,model,npts_per_fit)
	%matlabpool_cluster(4)
	debugmode = false;
	%npts_per_fit = 1e5;

	padding_minutes = 5; % Start this many minutes before sleep onset

	d = bt.load_seizure_data(subject_idx,'Cz',30,4); 
	st(6) = 1200;
	st(10) = 280;
	st(11) = 720;
	start_idx = st(subject_idx);
	d.t = d.t(start_idx:end);
	d.s = d.s(:,start_idx:end);
	d.s = bsxfun(@rdivide,d.s,trapz(d.f,d.s));

	% keyboard
	% for j = 1:1000
	% loglog(d.f,d.s(:,j))
	% pause
	% end

	model.set_electrodes('Cz');
	max_iters = 45; % 3600=1 hour
	initial_idx = 1;

	% figure
	% e_plot = loglog(d.f,d.s(:,1));
	% hold on
	% f_plot = loglog(d.f,d.s(:,1),'r--');
	% set(gca,'XLim',[1 45])
	
	for j = 1:max_iters
		tic;
		fprintf('Fitting: %d/%d\n',j,max_iters)
		
		target_P = squeeze(d.s(:,initial_idx+(j-1)));

		if j == 1
			[initial_params,initial_pp] = model.initialize_fit(d.f,target_P);
			[~,fit_data(j),plot_data(j)] = bt.core.fit_spectrum(model,d.f(:),target_P,initial_pp,initial_params,npts_per_fit,[],[],debugmode);
			fdata = bt.feather(model,fit_data(j),plot_data(j));
		else
			[~,fit_data(j),plot_data(j)] = bt.core.fit_spectrum(model,d.f(:),target_P,fit_data(j-1).posterior_pp,fit_data(j-1).fitted_params,npts_per_fit,[],[],debugmode);			
            fdata.insert(fit_data(j),plot_data(j));
		end

		%set(e_plot,'YData',target_P);
		%set(f_plot,'YData',fdata.fit_data(fdata.latest).fitted_P);
		%drawnow
		fdata.plot(j)
	end
	toc;
	fit = fdata;
	%save(sprintf('seizure_fit_test_%d.mat',subject_idx),'fit')

	% fdata.compress()
	% fdata.plot
