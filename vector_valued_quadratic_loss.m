function [L] = vector_valued_quadratic_loss(draws, central)
% Computes the quadratic loss for vector-valued draws relative to a central
% tendency measure.

ndraws = size(draws,3);
collapsed_draws = permute(draws, [3, 1, 2]);
collapsed_draws = collapsed_draws(:,:);

switch central
    case 'Bayes'
        L = NaN(ndraws,1);
        for i=1:ndraws
            L(i) = sum(sum((collapsed_draws-collapsed_draws(i,:)).^2));
        end
        L = L/ndraws;
    case 'mean'
        L = sum((collapsed_draws - mean(collapsed_draws)).^2,2);
    case 'median'
        L = sum((collapsed_draws - median(collapsed_draws)).^2,2);
end

end