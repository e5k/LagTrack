function download_SRTM(varargin)
% Function DOWNLOAD_SRTM - downloads tiles from the 90-m SRTM DEM dataset.
% Input either 0 or 6 arguments:
% lat_min:  Minimum latitude (decimal degree, negative in S hemisphere)
% lat_max:  Maximum latitude (decimal degree, negative in S hemisphere)
% lon_min:  Minimum longitude (decimal degree, negative in W hemisphere)
% lon_max:  Maximum longitude (decimal degree, negative in W hemisphere)
% res    :  Resolution (m) (Leave 90 m for no interpolation)
% name   :  File name, saved in input/dem/

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'Minimum latitude (decimal degree, negative in S hemisphere)', 'Maximum latitude (decimal degree, negative in S hemisphere)', 'Minimum longitude (decimal degree, negative in W hemisphere)', 'Maximum longitude (decimal degree, negative in W hemisphere)', 'Resolution (m)', 'DEM name'}, 'Download SRTM DEM', 1, {'','','','','90',''});
    if isempty(answer)
        return
    end
    lat_min     = str2double(answer{1});
    lat_max     = str2double(answer{2});
    lon_min     = str2double(answer{3});   
    lon_max     = str2double(answer{4});
    res         = str2double(answer{5});
    filename    = answer{6};   
elseif nargin == 6   
    lat_min     = varargin{1};
    lat_max     = varargin{2};
    lon_min     = varargin{3};   
    lon_max     = varargin{4};
    res         = varargin{5};
    filename    = varargin{6};
else
    error('Wrong number of input arguments, should be either 0 or 6');
end

% Check if folder already exist
if exist(['input/dem/', filename], 'dir') == 7
    choice = questdlg('A folder with the same name already exists. Do you want to retrieve a previous attempt or replace it?', ...
	'DEM Name', ...
	'Retrieve', 'Replace', 'Cancel', 'Retrieve');
    % Handle response
    switch choice
        case 'Replace' % Create a new run
            rmdir(['input/dem/', filename], 's');   % Remove existing folder
            mkdir(['input/dem/', filename]);        % Create folder
            
            % Retrieve tiles names and indices
            tiles = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max);
            
            dem.lat_min = lat_min;
            dem.lat_max = lat_max;
            dem.lon_min = lon_min;
            dem.lon_max = lon_max;
            dem.res     = res;
            dem.tiles   = tiles;
            dem.type    = 'DEM';
            save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');
            
        case 'Retrieve'
            load(['input/dem/', filename, filesep, filename, '.mat']);
            
        case 'Cancel'
            return
    end
end

% Main 2 loops
disp('Accessing SRTM server, please wait...')
maindir = ['input/dem/', filename];
for i = 1:length(dem.tiles)
    outdir  = [maindir, filesep, dem.tiles{i}];    % Tmp directory
    fprintf('   Downloading SRTM tile %s (%d of %d)... \n', tiles{i}, i, length(dem.tiles))
    
    % Check if tile is already downloaded
    if exist(outdir, 'dir')==7 && exist([outdir, filesep, tiles{i}, '.asc'], 'file')
        continue
    end
    
    % Else download
    DL_check = 0;
    while DL_check == 0
        DL_check = 1;
        websave([maindir, filesep, tiles{i}, '.zip'], ['http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_v41/SRTM_Data_ArcASCII/', tiles{i}, '.zip']);
        try
            unzip([maindir, filesep, tiles{i}, '.zip'], outdir);
        catch ME
            if strcmp(ME.identifier, 'MATLAB:unzip:invalidZipFile')
                DL_check = 0;
            end
        end
    end   
end

process_SRTM(lat_min, lat_max, lon_min, lon_max, res, filename);

