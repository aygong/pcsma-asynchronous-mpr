function r = natural_geornd(p, m, n)
% Generate random natural numbers from a geometric distribution
% with probability parameter p
r = zeros(m, n);
for i = 1:m
    for j = 1:n
        while r(i, j) == 0
            r(i, j) = geornd(p);
        end
    end
end