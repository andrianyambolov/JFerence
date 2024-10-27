function [out] = total_relative_interval_error(x1, x0)

tmp0 = squeeze(x0(:,:,2,:) - x0(:,:,1,:));
tmp1 = squeeze(x1(:,:,2,:) - x1(:,:,1,:));
out = mean(sum(sum((tmp1 - tmp0) ./ tmp0, 1),2));

end