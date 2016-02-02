function spectrum(self)
	% Plot the spectrum
	if self.latest == 0
		return
	elseif self.latest > 1
		fit = self.subrange(self.latest);
		fprintf(2,'Multiple fits provided - displaying the last one\n')
	else
		fit = self;
	end

	figure

	set(gca,'XScale','log','YScale','log')
	hold(gca,'on');
	box(gca,'on');
	set(gca,'XLim',[1 45]);
	%set(gca,'YLim',[3e-4 1]);
	xtf = [1:2:10 15 20 30 45];
	set(gca,'XTick',xtf);
	spec_exp = plot(self.fit_data.target_f,self.fit_data.target_P,'b');
	spec_fit = plot(self.fit_data.target_f,self.fit_data.fitted_P,'r');