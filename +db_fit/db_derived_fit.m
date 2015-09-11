function [pp_f,P_fitted,idx] = db_derived_fit(exp_f,exp_P,pp_f,pp,quiet)
    % Fit a power spectrum to a parameter combination already in a database
    % Uses derived spectral features rather than actual frequencies

    % Arguments are
    % f - experimental frequencies
    % P - experimental power
    % pp_f - database frequencies
    % pp - matrix of all power spectra in database
    if nargin < 6 || isempty(constraint_weights)
        constraint_weights = [];
    end

    if nargin < 5 || isempty(quiet)
        quiet = 1;
    end
    
    
    % Get the analysis
    exp_analysis = get_spec_analysis(exp_f,exp_P);

    features = {'int_delta','int_theta','int_alpha','int_sigma','int_beta','alpha_str','sigma_str','beta_str'};
    weights = [1 1 1 1 1 1 1 1]; 
    this_features = zeros(1,length(features));
    for j = 1:length(features)
        this_features(j) = exp_analysis.(features{j});
    end

    
    % Find the target frequency for normalization if applicable (currently unused)
    %pkf = get_pkf(pp_f,P);
    %target_index = find(pp_f == pkf(16)); % Target for normalization (spindle frequency)
    %target_index = find(pp_f == pkf(11)); % Target for normalization (alpha frequency)

    lsq = calc_derived_lsq(P,pp,weights);

    if ~isempty(constraint_weights)
        % Convert the LSQ into a probability
        % Then choose the most 'likely' fit
        prob = 1./lsq(:).*constraint_weights(:);
        [~,idx] = max(prob);
        % Next, apply a weighting based on XYZ
        
    else
        [lsq,idx] = sort(lsq);
    end

    P_fitted = pp(idx(1),:);
    normconst = P_fitted.' \ P.';
    P_fitted = P_fitted * normconst;

    %P_fitted = P_fitted * (P(target_index)/P_fitted(target_index));
    %P_fitted
    %keyboard

function lsq = calc_derived_lsq(P,pp,weights)
    sqdiff = bsxfun(@minus,pp,P).^2;
    lsq = sum(bsxfun(@times,sqdiff,weights),2);



