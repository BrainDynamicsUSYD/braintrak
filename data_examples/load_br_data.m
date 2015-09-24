function fdata = load_br_data(subject_id,electrode)
	if nargin < 2 || isempty(electrode)
		electrode = {'Cz'};
	elseif ischar(electrode)
		electrode = {electrode};
	end

	if nargin < 1 || isempty(subject_id)
		subject_id = '10002687';
	end

	if strcmp(subject_id,'demo')
		fdata = load('./braintrak/demo_spatial.mat');
	else
		fdata = load(sprintf('./psg_data/br_ec_multielectrode/%s.mat',subject_id));
	end

	% Find the indexes omatf fdata.electrodes that are present in electrode
	if ~isempty(electrode) && (length(electrode) > 1 || ~strcmp(electrode,'all'))
		idx = cellfun(@(x) find(strcmpi(x,fdata.electrodes)),electrode);
		fdata.electrodes = fdata.electrodes(idx);
		fdata.P = fdata.P(:,idx);
	end

	fidx = fdata.f >= 0.25 & fdata.f <= 45;
	fdata.f = fdata.f(fidx);
	fdata.P = fdata.P(fidx,:);

