function [coverage_error] = lqo_coverage_error(joint_cred, target_cred, draws, loss_value_idx, ndraws)
% Computes coverage error of loss quantile optimization.

cred_idxs = loss_value_idx(1:ceil(joint_cred*ndraws));
credible_set = draws(:,:,cred_idxs);
error_band = cat(3, min(credible_set, [], 3), max(credible_set, [], 3));
coverage_error = joint_coverage(draws, error_band) - target_cred;

end