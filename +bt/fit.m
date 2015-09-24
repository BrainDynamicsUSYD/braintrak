function f_out = fit(model,f,P,npts)
	% Minimal usage example 
	if nargin < 4 || isempty(npts)
		npts = [];
	end

	f_out =  bt.core.fit_spectrum(model,f,P,[],[],npts);
	f_out.plot()
