function [out] = relative_interval_area_difference(error_band1, error_band0)
% Computes the relative interval area difference using a harmonic mean.

nvars = size(error_band1,2);
area1 = NaN(nvars,1);
area0 = NaN(nvars,1);
ratio  = [];

for vv = 1:nvars
    tmp1 = abs(error_band1(:,vv,2) - error_band1(:,vv,1));
    tmp0 = abs(error_band0(:,vv,2) - error_band0(:,vv,1));
    ratio = [ratio; tmp1./tmp0];
end 

out = 100*(harmmean((ratio)) - 1);

end