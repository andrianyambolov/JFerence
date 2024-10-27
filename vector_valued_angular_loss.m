function [L] = vector_valued_angular_loss(draws, central)
% Computes the angular loss for vector-valued draws relative to a central
% tendency measure.

[~, nvars, ndraws] = size(draws);

if strcmp(central, 'Bayes')
        L = NaN(ndraws, nvars);
        for v = 1:nvars
            tmp = squeeze(draws(:,v,:));
            norm_collapsed_draws = (tmp ./ sqrt(sum(tmp.^2,1)))';
            Lvar = NaN(ndraws,ndraws);
            for i = 1:ndraws
                Lvar(:,i) = (sum(norm_collapsed_draws .* norm_collapsed_draws(i,:),2));
            end
            L(:,v) = sum(acos(Lvar),2);
        end
        L = sum(L,2) / (ndraws*nvars*pi);
else
        L = NaN(ndraws, nvars);
        for v = 1:nvars
            tmp = squeeze(draws(:,v,:));
            norm_collapsed_draws = (tmp ./ sqrt(sum(tmp.^2,1)))';
            if strcmp(central, 'mean')
                cent = mean(tmp,2)';
            elseif strcmp(central, 'median')
                cent = median(tmp, 2)';
            end
            norm_cent = cent / norm(cent);
            Lvar = sum(norm_collapsed_draws .* norm_cent,2);
            L(:,v) = acos(Lvar);
        end
        L = sum(L,2) / (ndraws*nvars*pi);
end

end