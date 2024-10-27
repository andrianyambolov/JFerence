function [pointwise] = pointwise_inference(post_draws, options)
% Computes pointwise credible intervals for posterior draws using specified
% credibility levels.

% Inputs:
%   - post_draws (double array of shape [horizons, nvars, ndraws]): 
%       Posterior draws of the parameters.
%       * horizons: Time dimension.
%       * nvars: Number of variables.
%       * ndraws: Number of posterior samples or draws.
%
%   - options (struct with named fields)
%      * Credibility: Credibility levels for credible sets, e.g., [0.68, 0.90].
%
% Outputs:
%   - pointwise (struct): Structure containing the results of pointwise inference.
%       * center: The median values of post_draws across draws for each variable and time horizon.
%       * error_band: Pointwise error bands at each specified credibility level.
%       * credibility: The list of credibility levels used.
%       * central_tendency: The central (median) tendency across draws.
%       * nvars, horizons, ndraws: Dimensions of posterior draws.

arguments
    post_draws(:,:,:) double
    options.Credibility (1,:) {mustBeFloat} = [0.68, 0.90];
end

credibility = options.Credibility;

[horizons, nvars, ndraws] = size(post_draws);

pointwise = struct;
pointwise.center = median(post_draws,3);

for cred = credibility
    cred_str = int2str(100*cred);
    quants = [(1-cred)/2, 1-(1-cred)/2];
    error_band = quantile(post_draws, quants, 3);
    pointwise.(['error_band', cred_str]) = error_band;
end

central_tendency = median(post_draws,3);

% set output values
pointwise.credibility = credibility;
pointwise.central_tendency = central_tendency;
pointwise.nvars = nvars;
pointwise.horizons = horizons;
pointwise.ndraws = ndraws;
end