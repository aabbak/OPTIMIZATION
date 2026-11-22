function main(algorithm)
% MAIN  Entry point: metaheuristic reference-point selection for hybrid-geoid
% determination. The optimiser only selects WHICH GNSS/levelling points are used;
% the estimation model is the classical least-squares corrector surface.
%
% USAGE
%   main            % runs PSO (default)
%   main('gwo')     % 'pso' | 'gwo' | 'abc' | 'tsa'
%
% Results are written to results_<algorithm>/ in the current folder. Run this
% from the repository root so the analysis script finds the result folders.

if nargin < 1 || isempty(algorithm), algorithm = 'pso'; end
rng('shuffle');

% --- Resolve paths so data loads regardless of the current folder ---
here     = fileparts(mfilename('fullpath'));     % .../src
repoRoot = fileparts(here);
addpath(here);                                   % optimiser + model functions
dataDir  = fullfile(repoRoot,'data');

ro = pi/180;

%% --- Read the gravimetric geoid model grid ---
grid = load(fullfile(dataDir,'gravimetric_geoid_model.xyz'));   % [lat lon N_grav]
lat = grid(:,1)'; lon = grid(:,2)'; ngr = grid(:,3)';
step = lon(2)-lon(1); latMin = min(lat); lonMin = min(lon); nGrid = length(ngr);
latAxis=[]; lonAxis=[]; Ngrid=[];
for i = 1:nGrid
    f = round((lat(i)-latMin)/step + 1);
    l = round((lon(i)-lonMin)/step + 1);
    latAxis(f) = lat(i);  lonAxis(l) = lon(i);  Ngrid(f,l) = ngr(i);
end

%% --- Candidate (training) and fixed test points ---
train = load(fullfile(dataDir,'gnss_leveling_training.txt'));   % [lat lon N_geom]
test  = load(fullfile(dataDir,'gnss_leveling_test.txt'));
latTrain=train(:,1); lonTrain=train(:,2); geoTrain=train(:,3);
latTest =test(:,1);  lonTest =test(:,2);  geoTest =test(:,3);

%% --- Interpolate gravimetric N once at the candidate and test points ---
nGravTrain = gravimetricInterpolation(latTrain,lonTrain, latAxis,lonAxis,Ngrid, step,latMin,lonMin, ro);
nGravTest  = gravimetricInterpolation(latTest, lonTest,  latAxis,lonAxis,Ngrid, step,latMin,lonMin, ro);

%% --- Fixed-data struct for the objective function ---
GEO.latTrain=latTrain; GEO.lonTrain=lonTrain; GEO.geoTrain=geoTrain; GEO.nGravTrain=nGravTrain;
GEO.latTest =latTest;  GEO.lonTest =lonTest;  GEO.geoTest =geoTest;  GEO.nGravTest =nGravTest;
GEO.ro = ro;
GEO.nParam = 5;     % corrector surface: 4, 5 or 7 (5 recommended for this region)
GEO.wAcc   = 0.8;   % accuracy weight (derived from RMSE)
GEO.wSr    = 0.2;   % sparsity weight (fewer points)
GEO.Adesign_test = correctorDesignMatrixAtGrid(latTest,lonTest, latAxis,lonAxis, step,latMin,lonMin, GEO.nParam, ro);

%% --- Algorithm settings (shared) ---
D = size(train,1);   % candidate count = selection-vector length
popSize = 20;
nRun    = 30;
maxFes  = 5000;

%% --- Run ---
switch lower(algorithm)
    case 'pso', pso(nRun,popSize,D,maxFes,GEO);
    case 'gwo', gwo(nRun,popSize,D,maxFes,GEO);
    case 'abc', abc(nRun,popSize,D,maxFes,GEO);
    case 'tsa', tsa(nRun,popSize,D,maxFes,GEO);
    otherwise,  error('algorithm must be one of: pso, gwo, abc, tsa.');
end

fprintf('\n=== %s done. Results written to results_%s/ ===\n', upper(algorithm), lower(algorithm));
end
