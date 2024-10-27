function [L] = vector_valued_absolute_loss(draws, central)
% Computes the absolute loss for vector-valued draws relative to a central
% tendency measure.

ndraws = size(draws,3);

collapsed_draws = permute(draws, [3, 1, 2]);
collapsed_draws = collapsed_draws(:,:);

switch central
    case 'Bayes'
        L = NaN(ndraws,1);
        for i=1:ndraws
            L(i) = sum(sum(abs(collapsed_draws-collapsed_draws(i,:))));
        end
        L = L/ndraws;
    case 'mean'
        L = sum(abs(collapsed_draws - mean(collapsed_draws)),2);
    case 'median'
        L = sum(abs(collapsed_draws - median(collapsed_draws)),2);
end
         
end