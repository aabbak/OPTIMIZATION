function Adesign = correctorDesignMatrixAtGrid(latP, lonP, latAxis, lonAxis, step, latMin, lonMin, nParam, ro)
% CORRECTORDESIGNMATRIXATGRID  Corner-weighted corrector design matrix for the
% test points, faithful to the classical grid-hybridise-then-interpolate scheme.
%
% The hybrid geoid is built on the model grid (N_grav + a*x) and then
% IDW-interpolated at each test point. Since this is linear in x:
%   N_hybrid(test) = IDW(N_grav) + ( IDW-weighted mean of corner design rows ) * x
% The second term's matrix (Adesign) is precomputed here, so the grid never has
% to be rebuilt during optimisation and the result matches the classical scheme.
ns = numel(latP);
tmp = correctorDesignMatrix(latAxis(1), lonAxis(1), nParam, ro);
u   = size(tmp,2);
Adesign = zeros(ns,u);
for i = 1:ns
    f = floor((latP(i)-latMin)/step + 1);
    l = floor((lonP(i)-lonMin)/step + 1);
    latLo = latAxis(f);   latHi = latAxis(f+1);
    lonLo = lonAxis(l);   lonHi = lonAxis(l+1);
    k1 = 6371*sqrt(((latLo-latP(i))*ro)^2 + (cos(latLo*ro)*(lonLo-lonP(i))*ro)^2);
    k2 = 6371*sqrt(((latHi-latP(i))*ro)^2 + (cos(latHi*ro)*(lonLo-lonP(i))*ro)^2);
    k3 = 6371*sqrt(((latHi-latP(i))*ro)^2 + (cos(latHi*ro)*(lonHi-lonP(i))*ro)^2);
    k4 = 6371*sqrt(((latLo-latP(i))*ro)^2 + (cos(latLo*ro)*(lonHi-lonP(i))*ro)^2);
    w  = [1/k1; 1/k2; 1/k3; 1/k4];  w = w/sum(w);
    % four corners in the same order as the IDW weights
    Ac = correctorDesignMatrix([latLo;latHi;latHi;latLo], [lonLo;lonLo;lonHi;lonHi], nParam, ro);
    Adesign(i,:) = w' * Ac;
end
end
