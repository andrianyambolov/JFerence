function [error_band, credible_set, opt_val] = pointwise_quantile_optimization(cred, draws)
% Implements pointwise quantile optimization of the sup-t method.

coverage_error_simple = @(x) pqo_coverage_error(x, cred, draws);
opt_val = secant_search(coverage_error_simple, [cred, 1]);
quants = [(1-opt_val)/2, 1-(1-opt_val)/2];
error_band = quantile(draws, quants, 3);
credible_set = enclosed_draws(draws, error_band);

end