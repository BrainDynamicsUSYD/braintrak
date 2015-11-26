function import_coma

	os_prefix = './psg_data/';

	for j = 1:1
		import(os_prefix,'coma',j);
	end

	for j = 1:1
		import(os_prefix,'control_coma',j);
	end

function import(os_prefix,data_set,k)
	infile = fullfile(os_prefix,data_set,sprintf('%s_%d.mat',data_set,k));
	fprintf('Loading: %s\n',infile);
	fhandle = load(infile);
	colheaders = fhandle.colheaders;
	for j = 1:length(colheaders)
		fprintf('%s_%d: %s electrode\n',data_set,k,colheaders{j});
		[t,f,s(:,:,j),nspec(:,j),n_reject(j)] = bt.data.get_tfs(fhandle.t(:),fhandle.data(:,j),30);
	end

	t_max_retain = floor(length(t)/30)*30;
	t = t(1:t_max_retain);
	t = t(:);
	s = s(:,1:t_max_retain,:);
	nspec = nspec(1:t_max_retain,:);

	state_str = cell(length(t),1);
	[state_str{:}] = deal(data_set);
	state_score = fhandle.epoch_score;

	outfile = fullfile(os_prefix,data_set,sprintf('%s_%d_tfs.mat',data_set,k));
	fprintf('Saving: %s\n',outfile);
	save(outfile,'t','f','s','state_score','state_str','nspec','n_reject','colheaders')
	fprintf('Done\n\n')