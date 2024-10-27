function [boundary_set, boundary_idx] = boundary_draws(draws, error_band)
% Selects the draws on the boundary.

collapsed_vals = permute(draws, [3, 1, 2]);
collapsed_vals = collapsed_vals(:,:);

lb = reshape(error_band(:,:,1), 1, []);
ub = reshape(error_band(:,:,2), 1, []);

[~, lb_idx] = mink(abs(collapsed_vals - lb), 1);
lb_idx = lb_idx(:);

[~, ub_idx] = mink(abs(collapsed_vals - ub), 1);
ub_idx = ub_idx(:);

boundary_idx = [lb_idx; ub_idx];
boundary_idx = unique(boundary_idx);
boundary_set = draws(:,:,boundary_idx);

end