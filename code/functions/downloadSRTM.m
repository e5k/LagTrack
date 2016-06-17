function [XX,YY,ZZ] = downloadSRTM(lat_min, lat_max, lon_min, lon_max, res)

if lon_min > 180
    lon_min = -360+lon_min;
    lon_max = -360+lon_max;
end



lon_minA = floor(lon_min);
lon_maxA = ceil(lon_max);
lat_minA = floor(lat_min);
lat_maxA = ceil(lat_max);

lat_grid = fliplr(-60:5:55);
lon_grid = -180:5:180;

[~, lat_minI] = min(abs(lat_grid-lat_minA));
[~, lat_maxI] = min(abs(lat_grid-lat_maxA));
[~, lon_minI] = min(abs(lon_grid-lon_minA));
[~, lon_maxI] = min(abs(lon_grid-lon_maxA));

lat_minI = 5;
lat_maxI = 5;
lon_minI = 39;
lon_maxI = 40;

% Retrieve data
nrows  = 6001;
ncols  = 6001;

XX = zeros(length(lat_maxI:lat_minI)*nrows, length(lon_minI:lon_maxI)*ncols); 
YY = zeros(size(XX));
ZZ = zeros(size(XX));

count  = 1;
county = 1;
idxY   = 1;
for yy = lat_maxI:lat_minI
    idxX    = 1;
    countx  = 1;
    for xx = lon_minI:lon_maxI
% 
%          tmpdir      = 'tmp_DEM_000000';
%          mkdir(tmpdir);
%          display(sprintf('Downloading file %d of %d... ', count, length(lat_maxI:lat_minI)*length(lon_minI:lon_maxI)))
%          
%          DL_check = 0;
%          while DL_check == 0
%              DL_check = 1;
%              websave('tmp.zip', ['http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_v41/SRTM_Data_ArcASCII/srtm_', num2str(xx, '%02d'), '_', num2str(yy, '%02d'), '.zip']);
%              try
%                 unzip('tmp.zip', tmpdir);
%              catch ME
%                  if strcmp(ME.identifier, 'MATLAB:unzip:invalidZipFile')
%                      DL_check = 0;
%                  end
%              end
%          end
%          [tmpX, tmpY, tmpZ] = readDEM([tmpdir, filesep, 'srtm_', num2str(xx, '%02d'), '_', num2str(yy, '%02d'), '.asc']);
         

        [tmpX, tmpY, tmpZ] = readDEM([ 'srtm_', num2str(xx, '%02d'), '_', num2str(yy, '%02d'), '.asc']);


         XX(idxY:idxY+nrows-1, idxX:idxX+ncols-1) = tmpX;
         YY(idxY:idxY+nrows-1, idxX:idxX+ncols-1) = tmpY;
         ZZ(idxY:idxY+nrows-1, idxX:idxX+ncols-1) = tmpZ;
         
%          rmdir([tmpdir,filesep],'s');
%          delete('tmp.zip');
%          count = count+ 1;
         
        idxX   = countx*ncols+1;
        countx = countx + 1;
    end
    idxY    = county*nrows+1;
    county  = county + 1;    
end

% Clean coordinates
[XX,YY] = meshgrid(linspace(XX(1,1),XX(1,end), size(XX,2)), linspace(YY(1,1),YY(end,1), size(XX,1)));


% Extract zone of interest
[~, lat_minI] = min(abs(YY(:,1)-lat_min));
[~, lat_maxI] = min(abs(YY(:,1)-lat_max));
[~, lon_minI] = min(abs(XX(1,:)-lon_min));
[~, lon_maxI] = min(abs(XX(1,:)-lon_max));


XX = XX(lat_maxI:lat_minI, lon_minI:lon_maxI);
YY = YY(lat_maxI:lat_minI, lon_minI:lon_maxI);
ZZ = ZZ(lat_maxI:lat_minI, lon_minI:lon_maxI);

% XX = XX(lat_minI:lat_maxI, lon_minI:lon_maxI);
% YY = YY(lat_minI:lat_maxI, lon_minI:lon_maxI);
% ZZ = ZZ(lat_minI:lat_maxI, lon_minI:lon_maxI);

% Interpolate
cellsize = XX(1,2) - XX(1,1);

[xq,yq] = meshgrid(XX(1,1):res*cellsize/90:XX(1,end),YY(end,1):res*cellsize/90:YY(1,1));
zq = interp2(XX,YY,ZZ,xq,yq);

dem.X   = xq; 
dem.Y   = yq; 
dem.Z   = zq;
dem.res = res;

save('dem.mat', 'dem');

