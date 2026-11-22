function fit = fitnessTransform(cost)
% FITNESSTRANSFORM  Standard ABC fitness mapping of a cost value to a positive
% fitness. Used by abc.m with cost = -objective so that a larger objective maps
% to a larger (always positive) fitness, keeping the selection probabilities valid.
fit = zeros(size(cost));
idx = cost >= 0;  fit(idx) = 1./(cost(idx)+1);
idx = cost <  0;  fit(idx) = 1 + abs(cost(idx));
end
