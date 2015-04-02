function output = load_seizure(subject_idx,electrode,window_length,fft_length) 
	% Generate a series of spectra from the seizure data
	% 6, 10, 11, 13 are good examples

	if nargin < 3 || isempty(electrode)
		electrode = {'Cz'};
	elseif ischar(electrode)
		electrode = {electrode};
	end

	if strcmp(electrode{1},'all')
		electrode = {'F4','Fp2','F3','Fp1','T6','T5','O2','O1','F7','F8','T3','T4','C4','C3','P4','P3','Cz','ECGp','Pz','A1','A2','Fz'};
	end

	if isunix && ~ismac
		fdata = load(sprintf('~/Desktop/psg_data/epilepsy/seizure_%d',subject_idx));
    elseif ismac
    	fdata = load(sprintf('~/Desktop/psg_data/epilepsy/seizure_%d',subject_idx));
    else
    	fdata = load(sprintf('C:/Users/romesh/Desktop/psg_data/epilepsy/seizure_%d',subject_idx));
    end

    t = (0:(size(fdata.raw_data,1)-1))/200;
    data = array2table([t(:) fdata.raw_data],'VariableNames',['Time',fdata.electrodes]);
    ref = (data.A1 + data.A2)/2;
    data{:,electrode }= bsxfun(@minus,data{:,electrode},ref); % re-reference
    
    for j = 1:length(electrode)
    	[output.t,output.f,output.s(:,:,j), output.nspec(:,j), output.n_reject(j)] = tracking.get_tfs(data.Time(:),data.(electrode{j}),window_length,fft_length,false);
    	%figure
    	%tracking.spectrogram_eeg_imagesc(output.t,output.f,output.s(:,:,j),data.Time,data.(electrode{j})) 
    end
    % figure
    % tracking.spectrogram_eeg_imagesc(fdata.t,fdata.f,fdata.s,data.Time,data.Cz) 

    % keyboard
    % tracking.animate_spectrum(output)
    % subplot(3,1,1)
    % title(sprintf('Seizure Subject %d',subject_idx));
    % pfig(sprintf('seizure_%d',subject_idx))
    % close all
    % output.colheaders = electrode;
    % size(output.t)
    % save seizure_test output

	tic;
