function [coverage_error] = pqo_coverage_error(pointwise_cred, target_cred, draws)
% Computes the coverage error of a pointwise quantile optimization method.

quants = [(1-pointwise_cred)/2, 1-(1-pointwise_cred)/2];
error_band = quantile(draws, quants, 3);
coverage_error = joint_coverage(draws, error_band) - target_cred;

end