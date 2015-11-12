function [f_fitted,P_fitted,idx,db_data,min_chisq] = quick_fit(model,target_f,target_P,force_wake_state)
    % Fit a power spectrum to a parameter combination already in a database
    % Arguments are
    % model - BrainTrak model
    % target_f - experimental frequencies
    % target_P - experimental power
    % force_wake_state - set to 1 for wake and 0 for sleep, to fit to only these states. Leave empty otherwise
    import utils.mex_trapz
    
    load db_fit/pp_allstates

    if nargin >= 4 && ~isempty(force_wake_state)
        filt = db_data.iswake == force_wake_state;
        db_data.P = db_data.P(filt,:);
        db_data.iswake = db_data.iswake(filt,:);
        db_data.xyz = db_data.xyz(filt,:);
        db_data.nus = db_data.nus(filt,:);
        db_data.gab = db_data.gab(filt,:);
    end

    target_P = (interp1(target_f,target_P,db_data.f,'linear','extrap'));
    weights = model.get_weights(db_data.f);
    normalization_target = mex_trapz(db_data.f(weights>0),target_P(weights>0));   

    chisq = zeros(1,size(db_data.P,1));
    norm_pp = db_data.P;
    for j = 1:size(db_data.P,1)
        norm_pp(j,:) = normalization_target*norm_pp(j,:)./mex_trapz(db_data.f(weights>0),norm_pp(j,weights>0));
        sqdiff = (abs(norm_pp(j,:)-target_P)./target_P).^2; % This is the squared fractional difference
        chisq(j) = sum(sqdiff(:).*weights(:));
    end

    [min_chisq,idx] = min(chisq); % min chisq is the same as max likelihood with no/uniform priors
    idx = idx(1);

    P_fitted = norm_pp(idx,:);
    f_fitted = db_data.f;

