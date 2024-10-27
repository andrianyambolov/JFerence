function [error_band, credible_set, opt_val] = loss_quantile_optimization(cred, draws, loss_value_idx)
% Performs loss quantile optimization.

ndraws = size(draws, 3);
coverage_error_simple = @(x) lqo_coverage_error(x, cred, draws, loss_value_idx, ndraws);
opt_val = secant_search(coverage_error_simple, [0.01, cred]);
cred_idxs = loss_value_idx(1:ceil(opt_val*ndraws));
credible_set = draws(:,:,cred_idxs);
error_band = cat(3, min(credible_set, [], 3), max(credible_set, [], 3));

end