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
			time_filter{j+1} = (1+(rate*j:(j+fft_length)*rate-1));
			[~,~,P(:,j+1)] = utils.rfft(V(time_filter{j+1}),rate);

			v_std(j+1) = std(V(time_filter{j+1}));

		end

		if remove_artifacts
			% Check clipping voltage
			clean_clipping = true(size(time_filter));
			for j = 1:length(time_filter)
				if sum(abs(V(time_filter{j}))>clipping) > 9 % Empirically, 9 is a good balance. This should be justified properly at some point
					clean_clipping(j) = false;
				end
			end

			% Check power standard deviations
			run_threshold = rate/2; % at least this many consecutive points for it to be an artifact
			in_range = @(pow,std_lim) pow < (mean(pow)+std(pow)*std_lim);
			hf_filter = fv > 30 & fv < 45; % hf cutoff is 4.5 Hz
			hf_pow = arrayfun(@(j) utils.mex_trapz(fv(hf_filter),P(hf_filter,j)),1:size(P,2));
			clean_range = in_range(v_std,4.5) & in_range(hf_pow,3); 

			% Check for runs of flat signals
			repeated_index = find(diff(V)==0); % The indexes that where i+1 = i
			runs = zeros(length(repeated_index),2);
			runs(1,:) = [repeated_index(1) repeated_index(1)];
			pointer = 1;
			for j = 2:length(repeated_index)
				if repeated_index(j) == runs(pointer,2)+1 % If this follows a run
					runs(pointer,2) = repeated_index(j);
				else
					pointer = pointer + 1;
					runs(pointer,:) = [repeated_index(j) repeated_index(j)];
				end
			end
			runs = runs(1:pointer,:);
			runs = [runs,runs(:,2)-runs(:,1)];
			runs = runs(runs(:,3)>run_threshold,:);
			run_reject = false(size(V));
			for j = 1:size(runs,1)
				run_reject(runs(j,1):runs(j,2)) = true;
			end

			% figure
			% plot(1:length(V),V)
			% for j = 1:size(runs,1)
			% 	set(gca,'XLim',[runs(j,1)-rate*3 runs(j,1)+rate*7])
			% 	set(gca,'XTick',[runs(j,1)-rate*3:rate:runs(j,1)+rate*7],'XTick',[runs(j,1)-rate*3:rate:runs(j,1)+rate*7],'XTickLabel',[0:10])
			% 	hold on
			% 	asdf = plot(runs(j,1):runs(j,2),V(runs(j,1):runs(j,2)),'r')
			% 	pause
			% 	delete(asdf)
			% end

			% Now go over the 4s windows again and see if they contained any run artifacts

			clean_runs = true(size(time_filter));
			for j = 1:length(time_filter)
				if any(run_reject(time_filter{j}))
					clean_runs(j) = false;
				end
			end
			
			clean = clean_clipping & clean_range & clean_runs;
			
		else
			clean = true(size(delta_pow));
		end

		% Now do the calculations for the 30 second spectra
		spectra = nan(length(fv),size(P,2)-n_avg);
		tv = 1:(size(P,2)-n_avg);
		
		nreject = sum(~clean);
		nspec = zeros(size(tv));

		for j = tv
			time_filter = j:(j+n_avg-1);
			time_filter = time_filter(clean(time_filter));
			nspec(j) = length(time_filter);
			if nspec(j) > 0
				spectra(:,j) = sum(P(:,time_filter),2);
			end
		end
		
		% And finally select only some of the frequencies and renormalize
		ffilt = fv > 0 & fv < 45;
		fv = fv(ffilt);
		spectra = spectra(ffilt,:);
		spectra = bsxfun(@rdivide,spectra,trapz(fv,spectra));

