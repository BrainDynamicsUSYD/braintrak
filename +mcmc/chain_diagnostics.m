function chain_diagnostics(model,out,posterior_out,j_burnin,automode)
	if nargin < 5 || isempty(automode)
		automode = false;
	end

	npoints = size(out,1);

	% Check the chain properties
	%a = round(0.2*npoints):round(0.3*npoints);
	a = j_burnin:(j_burnin+round(0.5*npoints));
	b = (j_burnin+round(0.5*npoints)):round(npoints);
	c = 1:j_burnin;

	m1 = mean(out(a,:));
	m2 = mean(out(b,:));
	m3 = mean(out(c,:));

	s1 = var(out(a,:));
	s2 = var(out(b,:));
	s3 = var(out(c,:));

	pcdiff = 100*abs((m1-m2)./m1);
	T = (m1-m2)./(sqrt(s1/length(a) + s2/length(b)));
	for j = 1:length(model.param_names)
		fprintf('%d\t %s \t%.3f%% difference\n',j,model.param_names{j},pcdiff(j));
	end

	for j = 1:length(model.param_names)
		fprintf('%d\t %s \t%.3f Geweke Z-score\n',j,model.param_names{j},T(j));
	end
	
	figure
	for j = 1:(size(out,2))
		subplot(size(out,2)+1,1,j)
		l = plot(1:j_burnin,out(1:j_burnin,j),'r');
		hold on
		l = plot(j_burnin:size(out,1),out(j_burnin:size(out,1),j),'k');
		yl = get(gca,'YLim');
		set(gca,'XLim',[1 size(out,1)]);
		patch([a(1) a(end) a(end) a(1)],[yl(1) yl(1) yl(2) yl(2)],0.85*ones(1,3),'LineStyle','none')
		patch([b(1) b(end) b(end) b(1)],[yl(1) yl(1) yl(2) yl(2)],0.7*ones(1,3),'LineStyle','none')
		pname = model.param_names{j};
		if pname(1)=='G'
			pname = ['G_{',pname(2:end),'}'];
		end
		ylabel(pname);
		uistack(l,'top')
		set(gca,'YLim',yl)
	end

	subplot(size(out,2)+1,1,j+1)
	l = plot(1:size(out,1),posterior_out,'k');
	ylabel('Posterior density')
	set(gca,'XLim',[1 size(out,1)]);

	xlabel('Step')

	initial_posterior = model.make_posterior(out(a,:));
	final_posterior = model.make_posterior(out(b,:));
	burnin_posterior = model.make_posterior(out(c,:));

	figure
	npars = model.n_params;	
	n_rows = ceil(npars/2);
	hist_ax = []
	for j = 1:model.n_params
		hist_ax(end+1) = subplot(n_rows,2,j);
		box(hist_ax(j),'on');
		hold(hist_ax(j),'on');
		%set(hist_ax(j),'YLim',[0 1],'XLim',[model.limits(1,j) model.limits(2,j)]);

		normconst = 1;%max(initial_posterior.y(:,j));
		plot(hist_ax(j),initial_posterior.x(:,j),initial_posterior.y(:,j)/normconst,'k--');
		normconst = 1;%max(final_posterior.y(:,j));
		plot(hist_ax(j),final_posterior.x(:,j),final_posterior.y(:,j)/normconst,'k');

		normconst = 1;%max(burnin_posterior.y(:,j));
		%plot(hist_ax(j),burnin_posterior.x(:,j),burnin_posterior.y(:,j)/normconst,'r:');

		hold(hist_ax(j),'off');
		set(hist_ax(j),'XLim',[model.limits(1,j) model.limits(2,j)]);
		axis tight
		xl = get(gca,'XLim');
		if xl(1) > 0 && xl(1)<1e-10
			xl(1) = 0;
		end
		if xl(2)<0 && xl(2)>-1e-10
			xl(2) = 0;
		end
		set(gca,'XLim',xl);
		pname = model.param_names{j};
		if pname(1)=='G'
			pname = ['G_{',pname(2:end),'}'];
		end
		ylabel(hist_ax(j),pname)
	end

	sum(diff(out(:,1))~=0)/size(out,1)

	if ~automode
		keyboard
	end