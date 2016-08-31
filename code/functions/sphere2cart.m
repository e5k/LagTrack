function [u,v,w] = sphere2cart(bearing, inclination, velocity)
% Bearing:      Degree from north
% Inclination:  Degree from vertical
% Velocity:     Initial velocity (m/s)

if bearing<0 || bearing>360
    error('Bearing should be comprised between 0 and 360');
elseif inclination<0 || inclination>90
    error('Inclination should be comprised between 0 and 90');
elseif velocity<0
    error('Velocity should be positive')
end
    
if bearing >90 && bearing < 360
    azimuth = 360-bearing+90;
elseif bearing>=0 && bearing <= 90
    azimuth = bearing*(-1)+90;
end

azimuth   = deg2rad(azimuth);
elevation = deg2rad(90-inclination);

[u,v,w] = sph2cart(azimuth, elevation, velocity);