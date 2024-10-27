function [coverage] = pointwise_coverage(vals, error_band)
% Computes pointwise coverage of error bands.

[horizons, nvars, ndraws] = size(vals);
nsim = size(error_band,4);

collapsed_vals = permute(vals, [3, 1, 2]);
collapsed_vals = collapsed_vals(:,:);
non_nan = ~isnan(collapsed_vals);

lb = squeeze3rd(error_band(:,:,1,:));
lb = permute(lb, [3, 1, 2]);
lb = lb(:,:);

ub = squeeze3rd(error_band(:,:,2,:));
ub = permute(ub, [3, 1, 2]);
ub = ub(:,:);

for cc = 1:size(collapsed_vals,2)
    coverage(1,cc) = mean(lb(non_nan(:,cc),cc) <= collapsed_vals(non_nan(:,cc),cc) & collapsed_vals(non_nan(:,cc),cc) <= ub(non_nan(:,cc),cc));
end

coverage = reshape(coverage, horizons, nvars);


end