function plot_timecourse(f,nocolor)
	if nargin < 2 || isempty(nocolor)
		nocolor = false;
	end

	m = f.model;
	if isempty(m.param_symbols)
		disp('Retrieving model parameter symbols and units using new instance of the class')
		m = feval(class(m));
	end

	param_names = m.param_names;
	xyz_req = ~m.xyz_present();
	xyz_names = {'X','Y','Z'};
	xyz_symbols = {'X','Y','Z'};
	xyz_units = {'','',''};
	xyz = f.xyz;
	
	yv = f.fitted_params;
	yv = [yv, xyz(:,xyz_req)];
	param_names = [param_names,xyz_names(xyz_req)];
	param_symbols = [m.param_symbols,xyz_symbols(xyz_req)];
	param_units = [m.param_units,xyz_units(xyz_req)];

	n_rows = ceil(length(param_names)/2);
	n_cols = 2;
	a_idx = @(a,b) n_cols*(a-1) + b;
	ax = [];

	xv = 1:f.latest;
	
	contaminated = tracking.chisq_outliers(f.chisq);
	yv(contaminated,:) = NaN;
	
	for j = 1:length(param_names)
		ax(end+1) = subplot(n_rows,n_cols,1+a_idx(ceil(j/2),mod(j-1,2)));
		box(ax(j),'on');

		if nocolor
			plot(xv,yv(:,j))
			hold on
		else
			tracking.plot_statecolored(xv,yv(:,j),f);
		end

		if isempty(param_units{j})
			ylabel(ax(j),param_symbols{j})
		else
			ylabel(ax(j),sprintf('%s (%s)',param_symbols{j},param_units{j}))
		end

		%ylabel(ax(j),m.param_names{j})
		drawnow
	end



