function db_data = db_fit_make_pp(state,cutdown)
	% This function returns the db_data corresponding to a single state
	% db_data contains the power spectrum in the matrix pp
	% it also has all the corresponding parameter values required to reconstruct it
	% The n_d cutdown can be applied at this point
	if nargin < 2 || isempty(cutdown)
		cutdown = 1; % Cutdown by default
	end

	states = {'eo','ec','rem','n1','n2','n3','n2s'};
	selected_state_index = find(strcmp(state,states));

	switch state
		case {'eo','ec'}
			load data/pdb_wake xyz_final gab_final nus_final
			iswake = 1;
		otherwise
			load data/pdb_sleep xyz_final gab_final nus_final
			iswake = 0;
	end

	load pdb_specvalid alllim

	xyz_final = xyz_final(alllim.(state),:);
	gab_final = gab_final(alllim.(state),:);
	nus_final = nus_final(alllim.(state),:);

	% Next, do a cutdown
	if cutdown
		keep = nd_cutdown(xyz_final,0.025,1,[0 -1 0]);
		xyz = xyz_final(keep,:);
		nus = nus_final(keep,:);
		gab = gab_final(keep,:);
	else
		xyz = xyz_final;
		nus = nus_final;
		gab = gab_final;
	end

	p = params(iswake); 
	p.gab = gab(end,:); % Last gab
	[~,pp_f,P] = dispersion_calc(p);

	pp = zeros(size(gab,1),length(pp_f));

	for j = 1:size(pp,1)
		p = params(iswake); 
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
	db_data.state = selected_state_index*ones(size(pp,1),1);
	db_data.states = states;
