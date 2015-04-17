function f2 = subrange(self,idxs)
	% Return a feather that corresponds to a subset/extraction of the 
	% input feather
	if ischar(self.plot_data{idxs(1)})
		f2 = mcmc.feather(self.model,self.fit_data(idxs(1)),fullfile(self.path_prefix,self.plot_data{idxs(1)}),self.time(idxs(1)));
	else
		f2 = mcmc.feather(self.model,self.fit_data(idxs(1)),self.plot_data{idxs(1)},self.time(idxs(1)));
	end

	if length(idxs)>=2
		for j = 2:length(idxs)
			if ischar(self.plot_data{idxs(j)})
				f2.insert(self.fit_data(idxs(j)),fullfile(self.path_prefix,self.plot_data{idxs(j)}),self.time(idxs(j)));
			else
				f2.insert(self.fit_data(idxs(j)),self.plot_data{idxs(j)},self.time(idxs(j)));
			end
		end
	end

	f2.compress();