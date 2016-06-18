function dem = downloadSRTM(lat_min, lat_max, lon_min, lon_max, res, filename)

% Check if folder already exist
if exist(['input/dem/', filename], 'dir') == 7
    choice = questdlg('A folder with the same name already exists. Do you want to delete it?', ...
	'DEM Name', ...
	'No', 'Yes', 'Yes');
    % Handle response
    switch choice
        case 'Yes'
            rmdir(['input/dem/', filename], 's');
        case 'No'
            return
    end
end

% Make output folder
mkdir(['input/dem/', filename])

% If longitudes are expressed in degrees west up to 360, correct them in
% negtive values
if lon_min > 180
    lon_min = -360+lon_min;
    lon_max = -360+lon_max;
end

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

%% Download data
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
        maindir = ['input/dem/', filename];
        tile    = ['srtm_', num2str(xx, '%02d'), '_', num2str(yy, '%02d')];
        outdir  = [maindir, filesep, tile];    % Tmp directory
        
        mkdir(outdir);
        
        display(sprintf('Downloading SRTM tile %d of %d... ', count, length(lat_maxI:lat_minI)*length(lon_minI:lon_maxI)))
        
        % Download
        DL_check = 0;
        while DL_check == 0
            DL_check = 1;
            websave([maindir, filesep, tile, '.zip'], ['http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_v41/SRTM_Data_ArcASCII/', tile, '.zip']);
            try
                unzip([maindir, filesep, tile, '.zip'], outdir);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:unzip:invalidZipFile')
                    DL_check = 0;
                end
            end
        end
        
        % Clean
        delete([maindir, filesep, tile, '.zip']);
        
        % Open the file and save data
        % Here, first flipud so the highest latitude is on index 1 in y
        % axis
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

YY = flipud(YY);

% Extract zone of interest
[~, lat_minI]   = min(abs(YY(:,1)-lat_min));
[~, lat_maxI]   = min(abs(YY(:,1)-lat_max));
[~, lon_minI]   = min(abs(XX(1,:)-lon_min));
[~, lon_maxI]   = min(abs(XX(1,:)-lon_max));

XX              = XX(lat_minI:lat_maxI, lon_minI:lon_maxI);
YY              = YY(lat_minI:lat_maxI, lon_minI:lon_maxI);
ZZ              = flipud(ZZ(lat_minI:lat_maxI, lon_minI:lon_maxI));
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
save(['input/dem/', filename, filesep, filename, '.dem'], 'dem');

% Plot
figure;
plot([min(dem.X(1,:)), max(dem.X(1,:))], [min(dem.Y(:,1)), max(dem.Y(:,1))], '.k'); 
[lonVec, latVec, imag] = plot_google_map('Maptype', 'terrain'); hold on;                           % Plot google background
surface(dem.X, dem.Y, dem.Z./1000, prepare_google_map(dem, lonVec, latVec, imag)); 
daspect([1,1,10]);
shading flat
axis tight
xlabel('Longitude');
ylabel('Latitude');
zlabel('Altitude (km)');

