function Nout = gravimetricInterpolation(latP, lonP, latAxis, lonAxis, Ngrid, step, latMin, lonMin, ro)
% GRAVIMETRICINTERPOLATION  Inverse-distance-weighted (IDW) interpolation of the
% gravimetric geoid model (regular grid) at arbitrary points.
%
%   latP,lonP      : latitude/longitude of the interpolation points (vectors)
%   latAxis,lonAxis: grid axis values (latitude, longitude)
%   Ngrid          : grid gravimetric geoid heights, Ngrid(fRow,lCol)
%   step           : grid spacing (degrees)
%   latMin,lonMin  : grid minimum latitude/longitude
%   ro             : pi/180
%
% The four surrounding grid nodes are weighted by the inverse of their
% great-circle distance (small-angle form, R = 6371 km).
ns   = numel(latP);
Nout = zeros(ns,1);
for i = 1:ns
    f = floor((latP(i)-latMin)/step + 1);     % grid cell of the point
    l = floor((lonP(i)-lonMin)/step + 1);
    latLo = latAxis(f);   latHi = latAxis(f+1);
    lonLo = lonAxis(l);   lonHi = lonAxis(l+1);
    k1 = 6371*sqrt(((latLo-latP(i))*ro)^2 + (cos(latLo*ro)*(lonLo-lonP(i))*ro)^2);
    k2 = 6371*sqrt(((latHi-latP(i))*ro)^2 + (cos(latHi*ro)*(lonLo-lonP(i))*ro)^2);
    k3 = 6371*sqrt(((latHi-latP(i))*ro)^2 + (cos(latHi*ro)*(lonHi-lonP(i))*ro)^2);
    k4 = 6371*sqrt(((latLo-latP(i))*ro)^2 + (cos(latLo*ro)*(lonHi-lonP(i))*ro)^2);
    Nout(i) = (Ngrid(f,l)/k1 + Ngrid(f+1,l)/k2 + Ngrid(f+1,l+1)/k3 + Ngrid(f,l+1)/k4) ...
            / (1/k1 + 1/k2 + 1/k3 + 1/k4);
end
end
