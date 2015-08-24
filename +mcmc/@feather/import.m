function f = import(dirname)
	disp('Loading file')
	d = load_fit(dirname);
	chisq = chisq_analysis(d.f,d.target_P,d.fitted_P,d.skip_fit);

	j = 1;

	[a,b] = get_data(d,chisq,dirname,j);
	m = mcmc.model.full;
	m.set_electrodes('Cz');
	f = mcmc.feather(m,a,b);

	disp('Iterating')
	for j = 2:length(d.t)
		[a,b] = get_data(d,chisq,dirname,j);
		f.insert(a,b);
	end
	f.compress();


function fit_data = load_fit(dirname)
	%try
		fit_data = load(sprintf('%s/fit_output',dirname));
		stop_idx = length(fit_data.p);
	% catch
	% 	fit_data = load(sprintf('%s/cached_output',dirname));
	% 	stop_idx = fit_data.initial_j;
	% end

	fit_data.state_str = fit_data.state_str(1:stop_idx);
	fit_data.target_P = fit_data.target_P(:,1:stop_idx);
	fit_data.fitted_P = fit_data.fitted_P(:,1:stop_idx);
	fit_data.p = fit_data.p(1:stop_idx);

	if ~isfield(fit_data,'skip_fit')
		fit_data.skip_fit = false(stop_idx,8);
	end
	
	if ~isfield(fit_data,'fit_mode')
		fit_data.fit_mode = 'full';
	end
	
	fit_data.xyz = reshape([fit_data.p.xyz],3,[]).';

function [fit_data,plot_data_fname] = get_data(d,chisq,dirname,j)
	fit_data.target_f = d.f;
	fit_data.target_P = d.target_P(:,j);
	fit_data.skip_fit = d.skip_fit(j,:);
	fit_data.fitted_params = [d.p(j).gabcd d.p(j).alpha(1) d.p(j).beta(1) d.p(j).t0 d.p(j).emg_a];
	fit_data.fitted_posterior = NaN;
	fit_data.posterior_pp = d.posterior_pp(j);
	fit_data.fitted_chisq = chisq(j);
	fit_data.fitted_P = d.fitted_P(:,j);
	fit_data.xyz = d.xyz(j,:);
	fit_data.fit_time = d.time_per_fit;
	fit_data.state_str = d.state_str{j};

	plot_data_fname = sprintf('%s/plot_data_%i',dirname,j);

function chisq_out = chisq_analysis(f,tp,fp,skip_fit)
	% Since fitted P is provided, it is already normalized
	m = mcmc.model.full;
    w = m.get_weights(f);
    
	for j = 1:size(tp,2)
		P = fp(:,j);
		target_P = tp(:,j);

		sqdiff = (abs(P-target_P)./target_P).^2; % This is the squared fractional difference
		chisq_out(j) = sum(sqdiff(:).*w(:));

	end

	chisq_out(skip_fit(:,1)>=3) = NaN; % Remove chisq when bad data was present



