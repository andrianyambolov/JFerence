function [coverage] = joint_coverage(vals, error_band)
nsim = size(error_band,4);

collapsed_vals = permute(vals, [3, 1, 2]);
collapsed_vals = collapsed_vals(:,:);
non_nan = ~any(isnan(collapsed_vals),2);

if nsim == 1
    lb = reshape(error_band(:,:,1), 1, []);
    ub = reshape(error_band(:,:,2), 1, []);

    cov_idxs = prod(lb <= collapsed_vals  & collapsed_vals  <= ub,2);
    coverage = mean(cov_idxs);
else
    lb = squeeze3rd(error_band(:,:,1,:));
    lb = permute(lb, [3, 1, 2]);
    lb = lb(:,:);

    ub = squeeze3rd(error_band(:,:,2,:));
    ub = permute(ub, [3, 1, 2]);
    ub = ub(:,:);

    cov_idxs = prod(lb(non_nan,:) <= collapsed_vals(non_nan,:)  & collapsed_vals(non_nan,:)  <= ub(non_nan,:),2);
    coverage = mean(cov_idxs);
end

end