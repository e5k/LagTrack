function [LAT, LON] = dist2ll(lat,lon,dx,dy)
% lat = old lat
% lon = old lon
% dx = x increm
% dy = y increm
% LAT and LON = new coordinates

R = 6371*10^3;

d = sqrt(dx^2+dy^2);
b = mod(atan2d(dx,dy),360);

LAT = asind(sind(lat) * cosd(d/R) + cosd(lat) * sin(d/R) * cosd(b));
LON = lon + atan2d( sind(b) * sin(d/R) * cosd(lat), cos(d/R) - sind(lat) * sind(LAT));