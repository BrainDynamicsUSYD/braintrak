function fdata = load_subject_data(dataset,subject_idx,electrode)
	if nargin < 3 || isempty(electrode)
		if strcmp(dataset,'control_apnea')
			electrode = {'Cz'};
		else
			electrode = {'C4-A1'};
		end
	end

	start_idx = 1;

	if any(strcmp({'control_apnea','control_opioid'},dataset))
		fdata = load(sprintf('./psg_data/sleep_tfs/%s_%d_hf_tfs',dataset,subject_idx));
		% For these data sets, skip the first bit of wake
		padding_minutes = 5; % Start this many minutes before sleep onset
		start_idx = max(1,find(fdata.state_score>2,1,'first')-padding_minutes*60);
	elseif any(strcmp({'olivia'},dataset))
		fdata = load(sprintf('./psg_data/olivia/%s_%d',dataset,subject_idx));
	else
		error('Data set unknown');
	end

	% Select the target electrodes
    if ~isempty(electrode) && (length(electrode) > 1 || ~strcmp(electrode,'all'))
   		idx = cellfun(@(x) find(strcmp(x,fdata.colheaders)),electrode);
   		fdata.colheaders = fdata.colheaders(idx);
   		fdata.nspec = fdata.nspec(:,idx);
   		fdata.n_reject = fdata.n_reject(idx);
   		fdata.s = squeeze(fdata.s(:,:,idx));
   	end

	fdata.t = fdata.t(start_idx:end);
	fdata.nspec = fdata.nspec(start_idx:end,:);
	fdata.s = fdata.s(:,start_idx:end,:);

	for j = 1:length(fdata.colheaders)
		fdata.s(:,:,j) = bsxfun(@rdivide,fdata.s(:,:,j),trapz(fdata.f,fdata.s(:,:,j)));
	end

	fdata.state_str = fdata.state_str(start_idx:end);
	fdata.state_score = fdata.state_score(start_idx:end);
	fdata.start_idx = start_idx;

