function assemble_master_pp
	states = {'eo','ec','rem','n1','n2','n3','n2s'};
	cutdown = [1 1 1 1 1 1 1 1];
	db_data = {};

	for j = 1:length(states)
		db_data{j} = db_fit_make_pp(states{j},cutdown(j));
	end

	for j = 2:length(states)
		db_data{1}.P = [db_data{1}.P; db_data{j}.P];
		db_data{1}.xyz = [db_data{1}.xyz; db_data{j}.xyz];
		db_data{1}.nus = [db_data{1}.nus; db_data{j}.nus];
		db_data{1}.gab = [db_data{1}.gab; db_data{j}.gab];
		db_data{1}.state = [db_data{1}.state; db_data{j}.state];
	end

	db_data = db_data{1};
	db_data.type = 'spec';


	% Next, work out the state transition constraint
	trans = count_transitions;
	t2 = zeros(7);
	t2(2:end,2:end) = trans; % Added an extra row and column for EO
	% Now, transitions TO EO can happen from any state going to wake, with equal probability of EO and EC
	% Transitions FROM EO can only go to EC (W)
	t2(1,:) = t2(2,:); % First row means that if a state goes TO EC, it could also go TO EO
	a = round(t2(2,2)/4);
	t2(1:2,1:2) = a;
	t2(1,2) = t2(3,2);
	t2 = bsxfun(@rdivide,t2,sum(t2,1)); % Normalized to transitions from

	db_data.trans = t2;

	save fitting/pp_master db_data

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
	scatter3(db_data.xyz(:,1),db_data.xyz(:,2),db_data.xyz(:,3),30,cdata(db_data.state,:))
