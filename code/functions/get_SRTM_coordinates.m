function [tiles,lat_minI, lat_maxI, lon_minI, lon_maxI] = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max)

% If longitudes are expressed in degrees east up to 360, correct them in
% negtive values
if lon_min > 180; lon_min = -360+lon_min; end
if lon_max > 180; lon_max = -360+lon_max; end

% Sort out SRTM coordinates
% Rounds input coordinates
lon_minA = floor(lon_min);
lon_maxA = ceil(lon_max);
lat_minA = floor(lat_min);
lat_maxA = ceil(lat_max);

% Vector of central points of srtm data
lat_grid = fliplr(-55+2.5:5:60-2.5);
lon_grid = -180+2.5:5:180-2.5;

[~, lat_minI] = min(abs(lat_grid-lat_minA));
[~, lat_maxI] = min(abs(lat_grid-lat_maxA));
[~, lon_minI] = min(abs(lon_grid-lon_minA));
[~, lon_maxI] = min(abs(lon_grid-lon_maxA));

% Create the tile list
tiles = cell(length(lat_maxI:lat_minI)*length(lon_minI:lon_maxI),1);
cnt = 1; % Counter
for yy = lat_maxI:lat_minI
    for xx = lon_minI:lon_maxI
        tiles{cnt} = ['srtm_', num2str(xx, '%02d'), '_', num2str(yy, '%02d')];
        cnt = cnt+1;
    end
end