function processSRTM(varargin)
% processSRTM Extract the region of interests from a set of 90-m SRTM tiles and interpolates to the desired resolution. 
%   processSRTM
%       Opens the GUI to select the root folder for all tiles.
%   processSRTM(lat_min, lat_max, lon_min, lon_max, res, name)
%       Process the SRTM tiles covering the specified extent
%           lat_min:  Minimum latitude (decimal degree, negative in S hemisphere)
%           lat_max:  Maximum latitude (decimal degree, negative in S hemisphere)
%           lon_min:  Minimum longitude (decimal degree, negative in W hemisphere)
%           lon_max:  Maximum longitude (decimal degree, negative in W hemisphere)
%           res    :  Resolution (m) (Leave 90 m for no interpolation)
%           name   :  File name, saved in input/dem/
%
%   see also downloadSRTM, makeDefaultGrid.
% 
% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'Minimum latitude (decimal degree, negative in S hemisphere)', 'Maximum latitude (decimal degree, negative in S hemisphere)', 'Minimum longitude (decimal degree, negative in W hemisphere)', 'Maximum longitude (decimal degree, negative in W hemisphere)', 'Resolution (m)'}, 'Process SRTM DEM', 1, {'','','','','90'});
    if isempty(answer)
        return
    end
    lat_min     = str2double(answer{1});
    lat_max     = str2double(answer{2});
    lon_min     = str2double(answer{3});   
    lon_max     = str2double(answer{4});
    res         = str2double(answer{5});
    
    % Select the folder
    folder = uigetdir('input/dem/', 'Select the folder containig the tiles');
    if isempty(folder)
        return
    end
    
    % Retrieve the file name
    tmp      = regexp(folder, filesep);
    filename = folder(tmp(end)+1:end);
    % If the .mat file is already present then load it, else create it
    if exist([folder, filesep, folder, '.mat'], 'file')
        load([folder, filesep, folder, '.mat']);
    else
        dem.lat_min = lat_min;
        dem.lat_max = lat_max;
        dem.lon_min = lon_min;
        dem.lon_max = lon_max;
        dem.res     = res;
        dem.tiles   = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max);
        dem.type    = 'DEM';
    end
    
elseif nargin == 6   
    lat_min     = varargin{1};
    lat_max     = varargin{2};
    lon_min     = varargin{3};   
    lon_max     = varargin{4};
    res         = varargin{5};
    filename    = varargin{6};
    % If the .mat file is already present then load it, else create it
    if exist(['input/dem/', filename, filesep, filename, '.mat'], 'file')
        load(['input/dem/', filename, filesep, filename, '.mat']);
    else
        dem.lat_min = lat_min;
        dem.lat_max = lat_max;
        dem.lon_min = lon_min;
        dem.lon_max = lon_max;
        dem.res     = res;
        dem.tiles   = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max);
        dem.type    = 'DEM';
        save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');
    end
else
    error('Wrong number of input arguments, should be either 0 or 6');
end


disp('Processing SRTM tiles, please wait...')

% Retrieves the indices of the files in the folder
maindir = ['input/dem/', filename];

% Retrieve indices of SRTM files
[~, lat_minI, lat_maxI, lon_minI, lon_maxI] = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max);

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
        [tmpX, tmpY, tmpZ] = readDEM([pwd, filesep, outdir, filesep, tile, '.hgt']);
        cellsize        = tmpX(1,2) - tmpX(1,1);            
        [xq,yq]         = meshgrid(tmpX(1,1):res*cellsize/90:tmpX(1,end),tmpY(1,1):res*cellsize/90:tmpY(end,1));
        zq              = interp2(tmpX,tmpY,tmpZ,xq,yq);
            
        % If first tile
        if idxX == 1 && idxY == 1
            % Define a storage matrix for the data of all tiles
            XX = zeros(length(lat_maxI:lat_minI)*size(xq,1), length(lon_minI:lon_maxI)*size(xq,1));
            YY = zeros(size(XX));
            ZZ = zeros(size(XX));
        end

        XX(idxY:idxY+size(xq,1)-1, idxX:idxX+size(xq,2)-1) = xq;
        YY(idxY:idxY+size(xq,1)-1, idxX:idxX+size(xq,2)-1) = flipud(yq);
        ZZ(idxY:idxY+size(xq,1)-1, idxX:idxX+size(xq,2)-1) = zq;

        count = count+ 1;
        
        % Update counters
        idxX   = countx*size(xq,2)+1;
        countx = countx + 1;
        clear tmpX tmpY tmpZ
    end
    idxY    = county*size(xq,1)+1;
    county  = county + 1;    
end

% Clean coordinates
[XX,YY] = meshgrid(linspace(XX(1,1),XX(1,end), size(XX,2)), linspace(YY(1,1),YY(end,1), size(XX,1)));
YY      = flipud(YY);
ZZ      = flipud(ZZ);

% Extract zone of interest
[~, lat_minI]   = min(abs(YY(:,1)-lat_min));
[~, lat_maxI]   = min(abs(YY(:,1)-lat_max));
[~, lon_minI]   = min(abs(XX(1,:)-lon_min));
[~, lon_maxI]   = min(abs(XX(1,:)-lon_max));

XX              = XX(lat_minI:lat_maxI, lon_minI:lon_maxI);
YY              = YY(lat_minI:lat_maxI, lon_minI:lon_maxI);
ZZ              = ZZ(lat_minI:lat_maxI, lon_minI:lon_maxI);
ZZ(isnan(ZZ))   = 0;

dem.X   = XX; 
dem.Y   = YY; 
dem.Z   = ZZ;
dem.res = res;

% Save
save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');

% Remove raw data
for yy = lat_maxI:lat_minI
    for xx = lon_minI:lon_maxI        
        tile    = ['srtm_', num2str(xx, '%02d'), '_', num2str(yy, '%02d')];
        if exist(['input/dem/', filename, filesep, tile], 'dir'); rmdir(['input/dem/', filename, filesep, tile], 's'); end
        if exist(['input/dem/', filename, filesep, tile, '.zip'], 'file'); delete(['input/dem/', filename, filesep, tile, '.zip']); end
    end
end

fprintf('Done!\n');