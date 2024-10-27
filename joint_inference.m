function [joint] = joint_inference(post_draws, options)
% Conducts joint inference about impulse responses and forecasts.
%
% Inputs:
%   - post_draws (double array of shape [horizons, nvars, ndraws]): 
%       Posterior draws of the parameters.
%       * horizons: Time dimension.
%       * nvars: Number of variables.
%       * ndraws: Number of posterior samples or draws.
%
%   - options: Specifies various options for inference.
%       * InferenceType (char): The type of inference to perform. 
%         Options are:
%           'MV' - Multivariate inference.
%           'TS' - Time-simultenous inference.
%           'MV-TS' - Combination of multivariate and time-simultenous inference.
%       * Method: The method for joint inference. 
%         Options are: 'min-max', 'sup-t', 'Bonferroni', 'Sidak'.
%       * LossFunction (optional): Specifies the loss function for 'min-max'.
%         Options are: 'absolute', 'quadratic', 'angular', 'chebyshev'.
%       * Calibration (optional): Calibration settings depending on the Method.
%       * Credibility: Credibility levels for credible sets, e.g., [0.68, 0.90].
%       * CentralTendency: Measure of central tendency. 
%         Options are: 'median', 'mean', 'Bayes'.
%
% Outputs:
%   - joint (struct): Contains results of the joint inference procedure.
%       * inference_type: Selected inference type.
%       * center: Computed central estimate (mean, median, or Bayesian).
%       * loss_value: Computed loss values if relevant.
%       * credible_set: Credible sets at different credibility levels.
%       * error_band: Computed error bands at different credibility levels.
%       * credibility: List of credibility levels used.
%       * loss_function: Specified loss function, if applicable.
%       * central_tendency: Measure of central tendency.
%       * is_calibrated: Indicates if calibration was applied.
%       * calibration: Calibration method, if applicable.
%       * nvars, horizons, ndraws: Dimensions of posterior draws.

arguments
    post_draws(:,:,:) double
    options.InferenceType (1,:) char {mustBeMember(options.InferenceType,{'MV','TS','MV-TS'})} = 'MV-TS';
    options.Method (1,:) char {mustBeMember(options.Method,{'min-max','sup-t','Bonferroni', 'Sidak'})} = 'min-max';
    options.LossFunction (1,:) char {mustBeText} = '';
    options.Calibration {mustBeText} = '';
    options.Credibility (1,:) {mustBeFloat} = [0.68, 0.90];
    options.CentralTendency (1,:) char {mustBeMember(options.CentralTendency,{'median','mean','Bayes'})} = 'median';
end

inference_type = options.InferenceType;
inference_method = options.Method;
loss_function = options.LossFunction;
calibration = options.Calibration;
credibility = options.Credibility;
central = options.CentralTendency;

% checks and balances
if strcmp(inference_method, 'min-max') && isempty(loss_function)
    loss_function = 'absolute';
end
loss_validation(inference_method, loss_function)

if ~isempty(calibration)
    calibration_validation(inference_method, calibration);
    is_calibrated = 1;
else
    is_calibrated = 0;
end

central_validation(inference_method, central)

if is_calibrated
    switch calibration
        case 'PQO'
            calibration_parameter_name = 'calibrated_pointwise_quantile';
        case 'LQO'
            calibration_parameter_name = 'calibrated_loss_quantile';
        case 'BDR'
            calibration_parameter_name = 'rejected_draws';
    end
end

% set output structure
[horizons, nvars, ndraws] = size(post_draws);

% output structure
joint = struct;
joint.inference_type = inference_type;

switch central
    case 'Bayes'
        joint.center = [];
        center_store = {};
    case 'mean'
        joint.center = mean(post_draws,3);
    case 'median'
        joint.center = median(post_draws,3);
end

if strcmp(inference_method, 'min-max')
    switch inference_type
        case 'MV-TS'
            joint.loss_value = NaN(1,1,ndraws);
        case 'TS'
            joint.loss_value = NaN(1,nvars,ndraws);
        case 'MV'
            joint.loss_value = NaN(horizons,1,ndraws);
    end
end

for cred = credibility
    cred_str = int2str(100*cred);
    joint.(['credible_set', cred_str]) =  NaN(horizons, nvars, ndraws);
    joint.(['error_band', cred_str]) = NaN(horizons, nvars, 2);
    if is_calibrated
        switch inference_type
            case 'MV-TS'
                joint.([calibration_parameter_name, cred_str]) = NaN(1, 1);
            case 'TS'
                joint.([calibration_parameter_name, cred_str]) = NaN(1, nvars);
            case 'MV'
                joint.([calibration_parameter_name, cred_str]) = NaN(horizons, 1);
        end
    end
end

%% conduct inference
flag = 0; sel = 0;
while flag < 1
    switch inference_type
        case 'MV-TS'
            sel = 1:nvars;
            post_draws_sel = post_draws(:,sel,:);
        case 'TS'
            sel = sel + 1;
            post_draws_sel = post_draws(:,sel,:);
        case 'MV'
            sel = sel + 1;
            post_draws_sel = post_draws(sel,:,:);
    end
    k = size(post_draws_sel, 1) * size(post_draws_sel, 2);

    switch inference_method
        case 'min-max'
            % evaluate loss function
            switch loss_function
                case 'absolute'
                    loss_value = vector_valued_absolute_loss(post_draws_sel, central);
                case 'quadratic'
                    loss_value = vector_valued_quadratic_loss(post_draws_sel, central);
                case 'angular'
                    loss_value = vector_valued_angular_loss(post_draws_sel, central);
                case 'chebyshev'
                    loss_value = vector_valued_chebyshev_loss(post_draws_sel, central);
            end

            [a, loss_value_idx] = sort(loss_value);
            loss_value = 100*(loss_value-a(1))/a(1);

            if strcmp(central, 'Bayes')
                c = log(ndraws)/sqrt(ndraws); % constant to estimate central tendency
                if strcmp(inference_type, 'MV-TS')
                    joint.center = post_draws_sel(:,:,loss_value<c);
                else
                    center_store{sel} = post_draws_sel(:,:,loss_value<c);
                end
            end
    end

    % loop over different credibility levels
    for cred = credibility
        cred_str = int2str(100*cred);

        switch inference_method
            case 'min-max'
                cred_idxs = loss_value_idx(1:ceil(cred*ndraws));
                credible_set = post_draws_sel(:,:,cred_idxs);
                error_band = cat(3, min(credible_set, [], 3), max(credible_set, [], 3));
            case 'Bonferroni'
                quants = [(1-cred)/(2*k), 1-(1-cred)/(2*k)];
                error_band = quantile(post_draws_sel, quants, 3);
                credible_set = enclosed_draws(post_draws_sel, error_band);
            case 'Sidak'
                quants = [(1-cred^(1/k))/2, 1-(1-cred^(1/k))/2];
                error_band = quantile(post_draws_sel, quants, 3);
                credible_set = enclosed_draws(post_draws_sel, error_band);
        end

        % in case of calibrated estimation
        if is_calibrated
            switch inference_method
                case 'min-max'
                    switch calibration
                        case 'BDR'
                            [error_band, credible_set, opt_val] = boundary_draw_rejection(cred, post_draws_sel, credible_set, error_band);
                        case 'LQO'
                            [error_band, credible_set, opt_val] = loss_quantile_optimization(cred, post_draws_sel, loss_value_idx);
                    end
                case 'Bonferroni'
                    [error_band, credible_set, opt_val] = boundary_draw_rejection(cred, post_draws_sel, credible_set, error_band);
                case 'Sidak'
                    [error_band, credible_set, opt_val] = boundary_draw_rejection(cred, post_draws_sel, credible_set, error_band);
                case 'sup-t'
                    [error_band, credible_set, opt_val] = pointwise_quantile_optimization(cred, post_draws_sel);
            end
            switch inference_type
                case 'MV-TS'
                    joint.([calibration_parameter_name, cred_str])(1) = opt_val;
                case 'TS'
                    joint.([calibration_parameter_name, cred_str])(:, sel) = opt_val;
                case 'MV'
                    joint.([calibration_parameter_name, cred_str])(sel, 1) = opt_val;
            end
        end

        % store credible set and error bands
        switch inference_type
            case 'MV-TS'
                joint.(['credible_set', cred_str])(:,:,1:size(credible_set,3)) =  credible_set;
                joint.(['error_band', cred_str])(:,:,:) = error_band;
            case 'TS'
                joint.(['credible_set', cred_str])(:,sel,1:size(credible_set,3)) =  credible_set;
                joint.(['error_band', cred_str])(:,sel,:) = error_band;
            case 'MV'
                joint.(['credible_set', cred_str])(sel,:,1:size(credible_set,3)) =  credible_set;
                joint.(['error_band', cred_str])(sel,:,:) = error_band;
        end
    end

    % store loss values and Bayes central estimates in case of min-max estimator
    % or stop the while if it is done
    switch inference_type
        case 'MV-TS'
            if strcmp(inference_method, 'min-max')
                joint.loss_value(1,1,:) = loss_value;
            end
            flag = 1;
        case 'TS'
            if strcmp(inference_method, 'min-max')
                joint.loss_value(1,sel,:) = loss_value;
            end
            if sel == nvars
                % set center in case of Bayes central tendency
                if strcmp(central, 'Bayes')
                    size_Bayes_store = [];
                    for isel = 1:numel(center_store)
                        size_Bayes_store = [size_Bayes_store, size(center_store{isel},3)];
                    end
                    max_size_Bayes = max(size_Bayes_store);
                    joint.center = NaN(horizons,nvars,max_size_Bayes);
                    for isel = 1:numel(center_store)
                        joint.center(:,isel,1:size(center_store{isel},3)) = center_store{isel};
                    end
                end
                flag = 1;
            end
        case 'MV'
            if strcmp(inference_method, 'min-max')
                joint.loss_value(sel,1,:) = loss_value;
            end
            if sel == horizons
                % set center in case of Bayes central tendency
                if strcmp(central, 'Bayes')
                    size_Bayes_store = [];
                    for isel = 1:numel(center_store)
                        size_Bayes_store = [size_Bayes_store, size(center_store{isel},3)];
                    end
                    max_size_Bayes = max(size_Bayes_store);
                    joint.center = NaN(horizons,nvars,max_size_Bayes);
                    for isel = 1:numel(center_store)
                        joint.center(isel,:,1:size(center_store{isel},3)) = center_store{isel};
                    end
                end
                flag = 1;
            end
    end
end

% set common output values
joint.credibility = credibility;
joint.loss_function = loss_function;
joint.central_tendency = central;
joint.is_calibrated = is_calibrated;
if is_calibrated
    joint.calibration = calibration;
end
joint.nvars = nvars;
joint.horizons = horizons;
joint.ndraws = ndraws;

end

%% validation functions
function central_validation(inference_method, central)
if ~strcmp(inference_method, 'min-max') && strcmp(central, 'Bayes')
    error('Joint Bayesian estimate of the central tendency is supported only for min-max estimator.')
end
end

function calibration_validation(inference_method, calibration)
switch inference_method
    case 'min-max'
        if ~(strcmp(calibration, 'BDR') || strcmp(calibration, 'LQO'))
            error('Min-max error band supports only BDR and LQO calibration.');
        end
    case 'sup-t'
        if ~strcmp(calibration, 'PQO')
            error('Sup-t error band supports only PQO calibration.');
        end
    case 'Bonferroni'
        if ~strcmp(calibration, 'BDR')
            error('Bonferroni error band supports only BDR calibration.');
        end
    case 'Sidak'
        if ~strcmp(calibration, 'BDR')
            error('Šidák error band supports only BDR calibration.');
        end
end
end

function loss_validation(inference_method, loss_function)
if strcmp(inference_method, 'min-max')
    if ~any(ismember({'absolute', 'quadratic', 'angular', 'chebyshev'}, loss_function))
        error([loss_function, ' is not supported.']);
    end
elseif ~strcmp(inference_method, 'min-max')
    if ~isempty(loss_function)
        error([inference_method, ' does not support loss function specifications.']);
    end
end
end

