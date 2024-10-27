function [out] = squeeze3rd(input)
% Squeeze the 3rd dimension of 4-dimensional matrix.

z = size(input);
out = reshape(input,[z(1:2), z(4), 1]);
end