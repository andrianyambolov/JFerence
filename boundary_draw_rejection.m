function [error_band, credible_set, opt_val] = boundary_draw_rejection(target_cred, draws, credible_set, error_band)
% Implements boundary draw rejection.

opt_val = 0;
coverage_error = joint_coverage(draws, error_band) - target_cred; 
while coverage_error > 0.0001
    [~, boundary_set_idx] = boundary_draws(credible_set, error_band);
    for bdraw = boundary_set_idx'
        credible_set(:,:,bdraw) = NaN;
    end
    error_band = cat(3, min(credible_set,  [], 3), max(credible_set,  [], 3));
    opt_val = opt_val + size(boundary_set_idx,1);
    coverage_error = joint_coverage(draws, error_band) - target_cred;
end

credible_set = enclosed_draws(draws, error_band);

end