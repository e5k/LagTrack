function process_SRTM(lat_min, lat_max, lon_min, lon_max, res, filename)

disp('Processing SRTM tiles, please wait...')

% Retrieves the indices of the files in the folder
maindir = ['input/dem/', filename];
files   = dir([maindir, filesep, 'srtm*']);
storI   = zeros(length(files),2);

for i = 1:length(files)
    storI(i,1) = str2double(files(i).name(6:7));
    storI(i,2) = str2double(files(i).name(9:10));
end

% Retrieve indices of SRTM files
[lat_minI, lat_maxI, lon_minI, lon_maxI] = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max);

% Size of each SRTM tile
nrows  = 6001;
ncols  = 6001;

% Define a storage matrix for the data of all tiles
XX = zeros(length(lat_maxI:lat_minI)*nrows, length(lon_minI:lon_maxI)*ncols); 
YY = zeros(size(XX));
ZZ = zeros(size(XX));

% Initialize counters
count  = 1;
county = 1;
idxY   = 1;

% Main 2 loops
for yy = lat_maxI:lat_minI
    idxX    = 1;
    countx  = 1;
    for xx = lon_minI:lon_maxI
        tile    = ['srtm_', num2str(xx, '%02d'), '_', num2str(yy, '%02d')];
        outdir  = [maindir, filesep, tile];    % Tmp directory
        % Open the file and save data
        % Here, first flipud so the highest latitude is on index 1 in y
        % axis
        disp(['   Reading tile ', tile, '...'])
        [tmpX, tmpY, tmpZ] = readDEM([pwd, filesep, outdir, filesep, tile, '.asc']);
        XX(idxY:idxY+nrows-1, idxX:idxX+ncols-1) = tmpX;
        YY(idxY:idxY+nrows-1, idxX:idxX+ncols-1) = flipud(tmpY);
        ZZ(idxY:idxY+nrows-1, idxX:idxX+ncols-1) = tmpZ;

        count = count+ 1;
        
        % Update counters
        idxX   = countx*ncols+1;
        countx = countx + 1;
    end
    idxY    = county*nrows+1;
    county  = county + 1;    
end

%% Post processing
% Clean coordinates
[XX,YY] = meshgrid(linspace(XX(1,1),XX(1,end), size(XX,2)), linspace(YY(1,1),YY(end,1), size(XX,1)));
YY      = flipud(YY);
ZZ      = flipud(ZZ);

% Extract zone of interest
[~, lat_minI]   = min(abs(YY(:,1)-lat_min));
[~, lat_maxI]   = min(abs(YY(:,1)-lat_max));
[~, lon_minI]   = min(abs(XX(1,:)-lon_min));
[~, lon_maxI]   = min(abs(XX(1,:)-lon_max));

disp('Done!');

XX              = XX(lat_minI:lat_maxI, lon_minI:lon_maxI);
YY              = YY(lat_minI:lat_maxI, lon_minI:lon_maxI);
ZZ              = ZZ(lat_minI:lat_maxI, lon_minI:lon_maxI);
ZZ(isnan(ZZ))   = 0;

% XX = XX(lat_maxI:lat_minI, lon_minI:lon_maxI);
% YY = YY(lat_maxI:lat_minI, lon_minI:lon_maxI);
% ZZ = ZZ(lat_maxI:lat_minI, lon_minI:lon_maxI);

% Interpolate
cellsize        = XX(1,2) - XX(1,1);
[xq,yq]         = meshgrid(XX(1,1):res*cellsize/90:XX(1,end),YY(1,1):res*cellsize/90:YY(end,1));
zq              = interp2(XX,YY,ZZ,xq,yq);

dem.X   = xq; 
dem.Y   = yq; 
dem.Z   = zq;
dem.res = res;

% Save
save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');