function D = ll2dist(lat1, lon1,lat2, lon2)
% lat1 = latitude of point 1
% lon1 = longitude of point 1


radius=6371;

lat1 = lat1*pi/180;
lat2 = lat2*pi/180;
lon1 = lon1*pi/180;
lon2 = lon2*pi/180;

deltaLat= lat2-lat1;
deltaLon= lon2-lon1;
a       = sin((deltaLat)./2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon./2).^2;
c       = 2.*atan2(sqrt(a),sqrt(1-a));
D       = radius.*c;    %Haversine distance
