function chain_diagnostics_parallel(model,out,chisq_out,accept)


	keyboard
	n = size(out{1},1);
	m = length(out);


	for k = 1:size(out{j},2) % For each variable
		for j = 1:m
			subplot(size(out{1},2),1,k);
			h = autocorr(out{j}(:,k),1000);
			plot(h)
			hold all
		end
	end

	theta = NaN;

	for k = 1:size(out{j},2) % For each variable
		for j = 1:m
			theta(j) = mean(out{j}(:,k));
			s2(j) = var(out{j}(:,k));
		end

		W = mean(s2);
		thetabar = mean(theta);
		B = n/(m-1)*sum((theta - thetabar).^2);
		V = (1-1/n)*W + 1/n*B;
		R(k) = sqrt(V/W);
	end
