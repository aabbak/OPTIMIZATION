function [metrics, objective, sel] = objectiveFunction(binSel, GEO)
% OBJECTIVEFUNCTION  Evaluate one candidate subset of GNSS/levelling points.
%
% The metaheuristic only decides WHICH reference points are used. For a given
% binary selection vector, a parametric corrector-surface (hybrid geoid) model
% is fitted by least squares on the selected points, and the external accuracy
% is measured on an independent test set. Fewer points with higher accuracy are
% rewarded. (No ANN/ELM: the estimation model is the least-squares corrector
% surface used in classical hybrid-geoid determination.)
%
% INPUT
%   binSel : 0/1 selection vector (length = number of candidate training points)
%   GEO    : struct with all fixed data (built once in main.m)
%            .latTrain .lonTrain .geoTrain .nGravTrain   (candidate points)
%            .latTest  .lonTest  .geoTest  .nGravTest    (fixed test points)
%            .Adesign_test  (corner-weighted design matrix for the test points)
%            .nParam .ro .wAcc .wSr
%
% OUTPUT
%   metrics : struct .RMSE .MSE .MAE .R .R2 .m0 .extRMS   (test-set metrics)
%   objective : scalar objective value (MAXIMISED)
%   sel     : indices of the selected candidate points

nParam = GEO.nParam;  ro = GEO.ro;

nAll = length(binSel);
sel  = find(binSel == 1);
nSel = length(sel);
sparsity = (nAll - nSel) / nAll;          % SR: grows as fewer points are used (0-1)

% --- Constraint: for least-squares adjustment the number of selected points
% must exceed the number of parameters by at least one (redundancy). Without
% redundancy the normal equations cannot be adjusted, so such a solution is
% infeasible and is penalised. i.e. nSel >= nParam + 1 (>= 6 for a 5-parameter model).
if nSel < nParam + 1
    metrics = struct('RMSE',Inf,'MSE',Inf,'MAE',Inf,'R',0,'R2',0,'m0',Inf,'extRMS',Inf);
    objective = -1;
    return;
end

% --- Fit the corrector surface on the selected points ---
latS = GEO.latTrain(sel);
lonS = GEO.lonTrain(sel);
lMeas = GEO.geoTrain(sel) - GEO.nGravTrain(sel);     % measurement vector: N_geom - N_grav
Asel  = correctorDesignMatrix(latS, lonS, nParam, ro);
x     = pinv(Asel' * Asel) * (Asel' * lMeas);        % x = (A'A)^-1 A'l

% --- Internal accuracy (a posteriori std, for reference only) ---
v = Asel*x - lMeas;
if nSel > nParam
    m0 = sqrt((v'*v)/(nSel - nParam));
else
    m0 = 0;
end

% --- External accuracy: hybrid geoid = gravimetric + corrector, at test points ---
% Faithful to the grid-hybridise-then-interpolate pipeline: the hybrid value at a
% test point is the IDW average of the four surrounding grid nodes. Because this
% is linear in x, the corner-weighted design matrix (Adesign_test) is precomputed
% once (in main.m) and only multiplied by x here.
Nhyb = GEO.nGravTest + GEO.Adesign_test * x;         % hybrid geoid at test points
res  = GEO.geoTest - Nhyb;                           % geometric - hybrid

% --- Error metrics (original scale: metres) ---
MSE  = mean(res.^2);
RMSE = sqrt(MSE);                 % zero-referenced RMSE (includes bias)
MAE  = mean(abs(res));
nt   = numel(res);
% External RMS as in classical hybrid-geoid evaluation (1-parameter / mean bias
% removed, divided by n-1). Directly comparable to the baseline paper's accuracy.
if nt > 1
    extRMS = sqrt(sum((res - mean(res)).^2)/(nt-1));
else
    extRMS = 0;
end

% R and R2
yt = GEO.geoTest(:);  yp = Nhyb(:);
if std(yt)==0 || std(yp)==0
    R = 0;
else
    c = corrcoef(yt,yp);  R = c(1,2);
    if isnan(R), R = 0; end
end
SS_res = sum(res.^2);
SS_tot = sum((yt - mean(yt)).^2);
if SS_tot==0, R2 = 0; else, R2 = 1 - SS_res/SS_tot; end

metrics = struct('RMSE',RMSE,'MSE',MSE,'MAE',MAE,'R',R,'R2',R2,'m0',m0,'extRMS',extRMS);

% --- Objective (RMSE-centred, maximised) ---
acc = 1/(1+RMSE);                 % approaches 1 as RMSE -> 0
objective = GEO.wAcc*acc + GEO.wSr*sparsity;
end
