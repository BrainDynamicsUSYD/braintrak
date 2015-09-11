function speed_racer(model)
	% Here he comes, here comes speed racer...
	runlength = [1e6,1e5,1e4,5e3,2e3];
	fnames = {'1e6','1e5','1e4','5e3','2e3'};
	for k = 1:length(fnames)
		n_runs = 20;
		f(n_runs) = mcmc.feather;
		parfor j = 1:n_runs
			f(j) = mcmc.speed_test(model,runlength(k));
		end

		save(sprintf('speed_racer_%s',fnames{k}));
	end

	% for j = 1:model.n_params
	% 	subplot(ceil(model.n_params/2),2,j);
	% 	for k = 1:n_runs
	% 		plot(arrayfun(@(x) x.fitted_params(j),f(k).fit_data))
	% 		hold on
	% 	end
	% 	xlabel(model.param_names{j});
	% 	set(gca,'YLim',model.limits(:,j));
	% 	%scatter(j,arrayfun(@(x) x.fit_data.fitted_params(j),f))
	% 	%hold on
	% end
