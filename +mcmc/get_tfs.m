function [tv,fv,spectra,nspec,nreject] = get_tfs(t,V,window_length,fft_length,remove_artifacts)
		% Calculate t,f,s based on input time and voltage
		% nspec is the number of spectra that were averaged 
		if nargin < 5 || isempty(remove_artifacts)
			remove_artifacts = true;
		end

		if nargin < 4 || isempty(fft_length)
			fft_length = 4;
		end

		if nargin < 3 || isempty(window_length)
			window_length = 30;
		end

		if max(abs(V)) > 500 % Some of the opioid data is not clipped to anything sensible
			clipping = 9999; % Disable the clipping check, since the data is just wierd
		else
			clipping = max(abs(V)) - 3; % Clip if it is within 3uV of the maximum. 155 for control, 125 for some of the opioid
		end

		rate = round(1/(t(2)-t(1))); % The sampling frequency. window*rate is the number of points needed
		time_step = 1; % 1 second difference between successive spectra
		
		stoptime = floor(t(end))-fft_length; % The time in seconds of th efstart of the last window
		[fv,~,P] = utils.rfft(V(1:rate*fft_length),rate);
		P = zeros(length(fv),stoptime+1);
		n_avg = (window_length-fft_length)/time_step+1; % Number of FFTs contributing to each spectrum

		clean = true(size(0:time_step:stoptime)); % Spectra all start out clean
		v_std = zeros(size(clean));
		% Get the set of 4s overlapping windows
		for j = 0:time_step:stoptime
			time_filter = (1+(rate*j:(j+fft_length)*rate-1));
			[~,~,P(:,j+1)] = utils.rfft(V(time_filter),rate);

			if sum(abs(V(time_filter))>clipping) > 9 % Empirically, 9 is a good balance. This should be justified properly at some point
				clean(j+1) = false;
				keyboard
			end

			v_std(j+1) = std(V(time_filter));

		end

		if remove_artifacts

			in_range = @(pow,std_lim) pow < (mean(pow)+std(pow)*std_lim);

			hf_filter = fv > 30 & fv < 45; % hf cutoff is 4.5 Hz
			hf_pow = arrayfun(@(j) utils.mex_trapz(fv(hf_filter),P(hf_filter,j)),1:size(P,2));

			clean = clean & in_range(v_std,4.5) & in_range(hf_pow,3); 

		else
			clean = true(size(delta_pow));
		end


		% Now do the calculations for the 30 second spectra
		spectra = zeros(length(fv),size(P,2)-n_avg);
		tv = 1:(size(P,2)-n_avg);
		
		nreject = sum(~clean);
		nspec = zeros(size(tv));

		for j = tv
			time_filter = j:(j+n_avg-1);
			time_filter = time_filter(clean(time_filter));
			nspec(j) = length(time_filter);
			spectra(:,j) = sum(P(:,time_filter),2);
		end
		
		% And finally select only some of the frequencies and renormalize
		ffilt = fv > 0 & fv < 45;
		fv = fv(ffilt);
		spectra = spectra(ffilt,:);
		spectra = bsxfun(@rdivide,spectra,trapz(fv,spectra));
