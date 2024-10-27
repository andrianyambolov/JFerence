function [out] = error_band_area(error_band1, error_band0)
% Computes the difference between the areas of two error bands.

nvars = size(error_band1,2);
area1 = NaN(nvars,1);
area0 = NaN(nvars,1);

for vv = 1:nvars
    area1(vv) = trapz(abs(error_band1(:,vv,2) - error_band1(:,vv,1)));
end 

for vv = 1:nvars
    area0(vv) = trapz(abs(error_band0(:,vv,2) - error_band0(:,vv,1)));
end 

out = sum(area1 - area0);

end