function import_raw_eeg
	% This function loads the mat file corresponding to the EDF data i.e.
	% e.g. psg_data\psgs\control_apnea\control_apnea_9\control_apnea_9.mat
	% It then performs Fourier transforms and artifact rejection, and then 
	% writes the output to psg_data/sleep_tfs

	os_prefix = './psg_data/';

	data_set = 'unresponsive';
	idx = 1;

	for j = idx
		import(os_prefix,data_set,j);
	end

function import(os_prefix,data_set,k)
	infile = fullfile(os_prefix,'psgs',data_set,sprintf('%s_%d',data_set,k),sprintf('%s_%d.mat',data_set,k));
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
	state_score = nan(1,length(t));

	states_available = {'W','R','1','2','3','4'};
	state_color_index = [2 3 4 5 6 6]; % Indexes for cdata

	for j = 1:length(t)
		if mod(j-1,30) == 0 % This is on a 30s block
			state_str{j} = sprintf('%s',fhandle.epoch_score{(j-1)/30+1});
			dominant_state = fhandle.epoch_score{(j-1)/30+1};
        else
			state_str{j} = sprintf('%s-%s (%d/%d)',fhandle.epoch_score{floor((j-1)/30)+1},fhandle.epoch_score{ceil((j-1)/30)+1},30-mod(j-1,30),mod(j-1,30));
			if 30-mod(j-1,30) > mod(j-1,30) % The dominant state is whichever one contributed more (time-wise)
				dominant_state = fhandle.epoch_score{floor((j-1)/30)+1};
			else
				dominant_state = fhandle.epoch_score{ceil((j-1)/30)+1};
			end
		end
		state_score(j) = state_color_index(strcmp(dominant_state,states_available));
	end
	outfile = fullfile(os_prefix,'sleep_tfs',sprintf('%s_%d_hf_tfs.mat',data_set,k));
	fprintf('Saving: %s\n',outfile);
	save(outfile,'t','f','s','state_score','state_str','nspec','n_reject','colheaders')
	fprintf('Done\n\n')