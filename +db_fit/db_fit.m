function [f_fitted,P_fitted,idx,prob] = db_fit(exp_f,exp_P,db_data,constraint_data)
    % Fit a power spectrum to a parameter combination already in a database
    % Arguments are
    % f - experimental frequencies
    % P - experimental power
    % pp_f - database frequencies
    % pp - matrix of all power spectra in database
    % constraint_data : 1/lsq*constraint_weight is the probability. constraint_weight can be calculated externally
    if nargin < 4 || isempty(constraint_data)
        constraint_data = [];
    end

    weights = get_weight_matrix(strcmp(db_data.type,'spec'),db_data.f);

    if strcmp(db_data.type,'spec')
        P = prepare_spectrum(exp_f,exp_P,db_data.f); % evaluate exp_P at the same frequencies as pp_f
        %pkf = get_pkf(pp_f,P);
        %target_index = find(pp_f == pkf(16)); % Target for normalization (spindle frequency)
        %target_index = find(pp_f == pkf(11)); % Target for normalization (alpha frequency)
        target_index = [];
        prob = calc_prob(P,db_data.P,weights,target_index); % Probability is large for most likely point

    else
        prob = calc_derived_prob(P,db_data.P,weights); % Probability is large for most likely point
    end

    if ~isempty(constraint_data)
        %constraint_data
        prob = reweight_probability(prob,constraint_data);
    end

    [~,idx] = max(prob);

    P_fitted = db_data.P(idx(1),:);
    normconst = P_fitted(:)\P(:);
    P_fitted = P_fitted * normconst;
    %P_fitted = P_fitted * (P(target_index)/P_fitted(target_index));
    %P_fitted
    %keyboard
    

    idx = idx(1);
    f_fitted = db_data.f;

function weights = get_weight_matrix(spec_mode,pp_f)
    if spec_mode
        weights = ones(size(pp_f));
        weights = pp_f.^-1; % Decades equally weighted
        weights(pp_f>25) = 0;
        %weights([1:4 7:22 25:length(weights)]) = 0; % Filter
    else
        features = {'int_delta','int_theta','int_alpha','int_sigma','int_beta','alpha_str','sigma_str','beta_str'};
        weights = [1 1 1 1 1 1 1 1]; 
    end


function prob = reweight_probability(prob,constraint_data)
    prob = constraint_data.xyz.*prob;

function P = prepare_spectrum(exp_f,exp_P,pp_f)
    % Interpolate the experimental spectrum onto the pp frequencies
    P = (interp1(exp_f,exp_P,pp_f,'linear','extrap'));

function prob = calc_prob(P,pp,weights,target_index)
    % This function returns the least squares matrix 
    
    % normconsts = arrayfun(@(idx) pp(idx,:)'\P',1:size(pp,1));
    % sqdiff = bsxfun(@minus,bsxfun(@times,pp,normconsts(:)),P).^2;
    % lsq = sum(bsxfun(@times,sqdiff,weights),2);

    for j = 1:size(pp,1)
        % normconst = pp(j,:)'\P';
        % pp(j,:) = pp(j,:) * normconst;
        % sqdiff = abs(pp(j,:)- P).^2;
        % lsq(j) = sum(sqdiff.*weights);

        normconst = lscov(P(:),pp(j,:).',weights);
        sqdiff = abs(pp(j,:)/normconst-P).^2;
        lsq(j) = sum(sqdiff.*weights);
    end



    %prob = 1./lsq(:);
    prob = exp(-lsq(:)); % Instead use a decaying exponential so that probability is correctly limited to 1

function prob = calc_derived_prob(P,pp,weights)
    sqdiff = bsxfun(@minus,pp,P).^2;
    lsq = sum(bsxfun(@times,sqdiff,weights),2);
    prob = exp(-lsq(:));
