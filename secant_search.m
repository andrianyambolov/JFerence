function [root] = secant_search(funct, bounds)
% Root-finding using secant search.

x0 = bounds(1);
x1 = bounds(2);

count = 0;
err = abs(x1-x0);
while err > 0.01 && count < 2000
    count = count + 1;
    x2=(x0*funct(x1)-x1*funct(x0))/(funct(x1)-funct(x0));
    if x2 > max(bounds)
        x2 = 1;
    elseif x2 < min(bounds)
        x2 = min(bounds);
    end
    x0 = x1;
    x1 = x2;
    err = abs(x1-x0);
    root = x2;
end

end