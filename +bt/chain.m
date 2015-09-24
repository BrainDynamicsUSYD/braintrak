function [out,final_posterior,accept_ratio] = chain(model,initial_values,n_points,debugmode,timelimit,n_burnin)
	% Take in a spectrum, prior distribution, and initial parameters
	% Return the MCMC chain and corresponding chisq values	

	% Burn-in method. First get an initial set of points, n_initial. Then take 
	% 10% of the requested points as burn in. Then return n_points samples from the chain
	% If n_burnin is specified, this is the number of points in initialization_stage_2
	% that are sampled. That is, the initial greedy sampling is still carried out
	% the same way


	% Non-greedy start used in Rev 1578
	if nargin < 5 || isempty(timelimit)
		timelimit = false;
	end

	if nargin < 4 || isempty(debugmode)
		debugmode = false;
	end

	model.validate()
	tic;
	

	% Preallocate
	n_initial = 100;

	if nargin < 6 || isempty(n_burnin)
		n_burnin = n_initial+round(n_points*0.1); % Total burn-in length
	else
		n_burnin = n_burnin + n_initial;
	end
	
	out = zeros(n_burnin+n_points,model.n_params);
	final_posterior = zeros(n_burnin+n_points,1);

	% Initialize the chain
	x = zeros(2,model.n_params);
	x(1,:) = initial_values; % Initial guess for the gains - this needs to be a valid point!
	out(1,:) = x(1,:);
	[chisq,current_posterior] = model.probability(initial_values);
	if ~isfinite(chisq)
		error('Initial point must have valid parameters')
	end
	accept = 1; % The first step in the chain was accepted by definition
	final_posterior(1) = current_posterior;
	

	fixed_cov = diag(0.05*model.initial_step_size);
	lambda = (2.38)^2/model.n_fitted;
	lambda = lambda/4;
	
	target_accept = 0.238;
	alpha = NaN;
	mu = NaN;

	% See page 15-16 of "Optimal proposal distributions and adaptive MCMC by Rosenthal 2011" 

	% To avoid numerical errors, only use the incremental covariance after n_initial steps have been taken
	
	initialization_stage = 1; % 1 - fixed step. 2 - cov burn in
	j = 2; % Replicates for loop, at least at the beginning
	j_stop = n_points;
	burnin_counter = 1;
	stage_1_time = 0;

	while j <= n_burnin+n_points
		% For speedtest, limit the execution time to achieve real-time performance
		if timelimit && mod(j,1000) % Only poll for the time limit
			elapsed = toc;

			if elapsed > timelimit 
				break
			elseif initialization_stage == 2 && elapsed > timelimit*0.2
				n_burnin = j; % Terminate the burn-in early
			end
		end

		if initialization_stage == 1
			x(2,:) = x(1,:) + model.initial_step_size.*randn(1,model.n_params); 
		else
			x(2,:) = mvnrnd_fast(x(1,:),lambda*my_cov);
		end
			
		x(2,model.skip_fit==2) = initial_values(model.skip_fit==2);

		[chisq,new_posterior] = model.probability(x(2,:));

		if isfinite(chisq)			
			alpha = min(1,new_posterior/current_posterior);

			if rand(1) < alpha % Note rand(1) <1 so if alpha == 1 it is automatically accepted
				x(1,:) = x(2,:);
				accept = accept+1;
				current_posterior = new_posterior;
			end
		end
	
		final_posterior(j) = current_posterior;
		out(j,:) = x(1,:);

		if initialization_stage == 1 && accept == n_initial % Activate covariance matrix
			% accept==n_initial means that j was set to n_initial AND that the point was accepted on this round
			initialization_stage = 2;
			stage_1_time = toc;
			my_cov = cov(out(1:j,:));
			mu = mean(out(1:j,:));
		elseif initialization_stage == 2 && j == n_burnin % Reset acceptance ratio counter
			initialization_stage = 3;
			elapsed = toc;
			stage_2_time = elapsed-stage_1_time;

			if timelimit
				eta = timelimit-elapsed;
			else
				eta = stage_2_time/(n_burnin-n_initial)*n_points;
			end
			fprintf('Burn-in finished: S1 %.2fs (%d pts), S2 %.2fs, Remain: ~%.2fs, Total: ~%.2fs\n',stage_1_time,burnin_counter,stage_2_time,eta,eta+elapsed);
			accept = 0; % Reset acceptance counter. The final acceptence ratio is then accept/n_points
		end

		if initialization_stage >= 2 % Do adaptation
			gammai = 1/(j);
			z_res = x(1,:)-mu;
			mu = mu + gammai*(z_res); % This should equal mean(out(1:j-1,:))
			my_cov = (j-2)/(j-1)*my_cov + gammai*(z_res'*z_res); % see wikipedia native online incremental variance formula
			lambda = exp( log(lambda) +  gammai*(alpha - target_accept));
			j = j + 1;
		else
			j = accept + 1; % Set j to be the next point - the greedy burnin is mentioned in Haario 2001
			burnin_counter = burnin_counter + 1;
			if burnin_counter == 1e8
				error('Initial burnin failed after 1000000 steps - initial step size must be incorrect')
			end
		end

	end

	if ~(initialization_stage==3)
		error('Burnin did not finish, something went wrong');
        keyboard
    end
	
	accept_ratio = accept/n_points;

	if debugmode
		bt.chain_diagnostics(model,out(1:(j-1),:),final_posterior(1:(j-1)),n_burnin);
	end
	
	out = out(n_burnin+1:(j-1),:);
	final_posterior = final_posterior(n_burnin+1:(j-1));

	fprintf('Accept: %.2f%% (%.2fs, %d points)\n',sum(accept/(j-1))*100,toc,(j-1)-n_burnin);

function [r,T] = mvnrnd_fast(mu,sigma,cases,T)
    [T,err] = cholcov_fast(sigma);
    r = randn(1,size(T,1)) * T + mu;

function [T,p] = cholcov_fast(Sigma)

    [T,p] = chol(Sigma);

    if p > 0
    	[n,m] = size(Sigma);

        [U,D] = eig(full((Sigma+Sigma')/2));

        % Pick eigenvector direction so max abs coordinate is positive
        [ignore,maxind] = max(abs(U),[],1);
        negloc = (U(maxind + (0:n:(m-1)*n)) < 0);
        U(:,negloc) = -U(:,negloc);

        D = diag(D);
        tol = eps(max(D)) * length(D);
        t = (abs(D) > tol);
        D = D(t);
        p = sum(D<0); % number of negative eigenvalues

        if (p==0)
            T = diag(sqrt(D)) * U(:,t)';
        else
            T = [];
            error('p~=0 after flag loop')
        end
    end

