function A = correctorDesignMatrix(lat, lon, nParam, ro)
% CORRECTORDESIGNMATRIX  Design matrix of the parametric corrector surface used
% in hybrid-geoid determination.
%   nParam = 4, 5 or 7
%   lat,lon : latitude/longitude in degrees, ro = pi/180
lat = lat(:);  lon = lon(:);
n = numel(lat);
switch nParam
    case 4   % 4-parameter model
        A = [cos(lat*ro).*cos(lon*ro), cos(lat*ro).*sin(lon*ro), sin(lat*ro), ones(n,1)];
    case 5   % 5-parameter model
        A = [cos(lat*ro).*cos(lon*ro), cos(lat*ro).*sin(lon*ro), sin(lat*ro), ...
             sin(lat*ro).*sin(lat*ro), ones(n,1)];
    case 7   % 7-parameter model
        e2 = 0.006694380023;                       % first eccentricity squared
        W  = sqrt(1 - e2*sin(lat*ro).^2);
        A = [cos(lat*ro).*cos(lon*ro), cos(lat*ro).*sin(lon*ro), sin(lat*ro), ...
             cos(lat*ro).*sin(lat*ro).*cos(lon*ro)./W, ...
             cos(lat*ro).*sin(lat*ro).*sin(lon*ro)./W, ...
             sin(lat*ro).*sin(lat*ro)./W, ones(n,1)];
    otherwise
        error('correctorDesignMatrix: nParam must be 4, 5 or 7.');
end
end
