function contaminated = chisq_outliers(chisq)
	if ~isnumeric(chisq) % If the user put a data object in instead
		chisq = chisq.chisq();
	end
	contaminated = utils.romesh_gesd(chisq,0.05,ceil(length(chisq)*0.1),1);
