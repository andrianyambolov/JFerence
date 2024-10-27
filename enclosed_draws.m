function [enclosed_set, enclosed_idx] = enclosed_draws(draws, error_band)
% Selects the draws contained within an error band.

collapsed_vals = permute(draws, [3, 1, 2]);
collapsed_vals = collapsed_vals(:,:);

lb = reshape(error_band(:,:,1), 1, []);
ub = reshape(error_band(:,:,2), 1, []);

enclosed_idx = find(prod(lb <= collapsed_vals  & collapsed_vals  <= ub, 2));
enclosed_set = draws(:,:,enclosed_idx);

end