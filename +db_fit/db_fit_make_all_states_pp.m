function varargout = db_fit_make_all_states_pp
	% This function returns the db_data corresponding to a single 
	% db_data contains the power spectrum in the matrix pp
	% it also has all the corresponding parameter values required to reconstruct it
	% The n_d cutdown can be applied at this point
	
	wake_data= load('data/pdb_wake');
	sleep_data= load('data/pdb_sleep');
	
	load data/pdb_specvalid alllim
	wake_valid = alllim.ec | alllim.eo;
	sleep_valid = alllim.rem | alllim.n1 | alllim.n2 | alllim.n3 | alllim.n2s;

	xyz_final = [wake_data.xyz_final(wake_valid,:); sleep_data.xyz_final(sleep_valid,:)];
	gab_final = [wake_data.gab_final(wake_valid,:); sleep_data.gab_final(sleep_valid,:)];
	nus_final = [wake_data.nus_final(wake_valid,:); sleep_data.nus_final(sleep_valid,:)];
	iswake = [ones(sum(wake_valid),1); zeros(sum(sleep_valid),1)];

	clear wake_data sleep_data

	% Now perform the nd cutdown
	% Note that if too few points are retained, there might not be enough
	% wake and sleep points in the regions where these states overlap
	keep = nd_cutdown(xyz_final,0.025,1,[0 -1 0]);
	xyz = xyz_final(keep,:);
	nus = nus_final(keep,:);
	gab = gab_final(keep,:);
	iswake = iswake(keep);

	p = params(iswake(end)); 
	p.gab = gab(end,:); % Last gab
	[~,pp_f,P] = dispersion_calc(p);

	pp = zeros(size(gab,1),length(pp_f));

	for j = 1:size(pp,1)
		p = params(iswake(j)); 
		p.gab = gab(j,:); % Last gab
		[~,~,pp(j,:)] = dispersion_calc(p);
	end

	% Exclude low frequencies to avoid spindle rolloff
	pp_f = pp_f(3:end);
	pp = pp(:,3:end);

	db_data.P = pp;
	db_data.f = pp_f;
	db_data.xyz = xyz;
	db_data.nus = nus;
	db_data.gab = gab;
	db_data.iswake = iswake;
	db_data.type = 'spec';

	save fitting/pp_allstates db_data

	if nargout > 0
		varargout{1} = db_data;
	end

	% Plot some output
	cdata = [255 0 0; ...% EO bright red
	0 255 145; ...% EC light green
	255 180 0; ...% REM Orange
	0 159 100; ...% N1 dark green
	0 154 255; ...% N2 light blue
	12 23 178; ...% N3 dark blue
	255 255 0; ...% N2S bright yellow
	];
	cdata = cdata/255;

	figure
	analysis_to_tent(0,[],0.1)
	scatter3(db_data.xyz(:,1),db_data.xyz(:,2),db_data.xyz(:,3),30)
