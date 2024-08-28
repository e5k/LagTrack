function downloadSRTM(varargin)
% downloadSRTM Download tiles from the 90-m SRTM DEM dataset (http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp) in ArcInfo ASCII format.
%   downloadSRTM
%       Open the GUI to specify the extent
%   downloadSRTM(lat_min, lat_max, lon_min, lon_max, res, name)
%       Download the SRTM tiles covering the specified extent
%           lat_min:  Minimum latitude (decimal degree, negative in S hemisphere)
%           lat_max:  Maximum latitude (decimal degree, negative in S hemisphere)
%           lon_min:  Minimum longitude (decimal degree, negative in W hemisphere)
%           lon_max:  Maximum longitude (decimal degree, negative in W hemisphere)
%           res    :  Downsampling factor at 1/res the resolution (1=30m, 3=90m)
%           name   :  File name, saved in input/dem/
%
%   see also processSRTM, makeDefaultGrid.
% 
% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'Minimum latitude (decimal degree, negative in S hemisphere)', 'Maximum latitude (decimal degree, negative in S hemisphere)', 'Minimum longitude (decimal degree, negative in W hemisphere)', 'Maximum longitude (decimal degree, negative in W hemisphere)', 'Downsampling factor (1/n)', 'DEM name'}, 'Download SRTM DEM', 1, {'','','','','1',''});
    if isempty(answer)
        return
    end
    lat_min     = str2double(answer{1});
    lat_max     = str2double(answer{2});
    lon_min     = str2double(answer{3});   
    lon_max     = str2double(answer{4});
    res         = str2double(answer{5});
    filename    = answer{6};
%     login       = answer{7};
%     passwd      = answer{8};
elseif nargin == 6   
    lat_min     = varargin{1};
    lat_max     = varargin{2};
    lon_min     = varargin{3};   
    lon_max     = varargin{4};
    res         = varargin{5};
    filename    = varargin{6};
%     login       = varargin{7};
%     passwd      = varargin{8};
else
    error('Wrong number of input arguments, should be either 0 or 6');
end

% Retrieve tiles names and indices 
%[tiles,~,~,~,~,LT,LN] = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max);

% Check if folder already exist
if exist(['input/dem/', filename], 'dir') == 7
    choice = questdlg('A folder with the same name already exists. Do you want to overwrite it?', ...
	'DEM Name', ...
	 'Replace', 'Cancel', 'Cancel');
    % Handle response
    switch choice
        case 'Replace' % Create a new run
            rmdir(['input/dem/', filename], 's');   % Remove existing folder
            mkdir(['input/dem/', filename]);        % Create folder

            dem.lat_min = lat_min;
            dem.lat_max = lat_max;
            dem.lon_min = lon_min;
            dem.lon_max = lon_max;
            %dem.res     = res;
            %dem.tiles   = tiles;
            %dem.lat     = LT;
            %dem.lon     = LN;
            dem.type    = 'DEM';
            save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');
            
%         case 'Retrieve'
%             load(['input/dem/', filename, filesep, filename, '.mat']);
%             
        case 'Cancel'
            return
    end
else
    mkdir(['input/dem/', filename]);        % Create folder    
    dem.lat_min = lat_min;
    dem.lat_max = lat_max;
    dem.lon_min = lon_min;
    dem.lon_max = lon_max;
    %dem.res     = res;
    %dem.tiles   = tiles;
    %dem.lat     = LT;
    %dem.lon     = LN;
    dem.type    = 'DEM';
    save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');    
end

% Check if folder _rawdata exists
if ~exist('input/dem/_rawdata', 'dir') == 7
    mkdir('input/dem/_rawdata');
end

% Download tiles
% demTmp  = readhgt([dem.lat_min, dem.lat_max, dem.lon_min, dem.lon_max], 'interp', 'outdir', fullfile(pwd,'input','dem','_rawdata'), 'srtm1', 'login', login, passwd);
demTmp  = readhgt([dem.lat_min,dem.lat_max, dem.lon_min, dem.lon_max], ...
    'interp', 'decim', res, ...
    'outdir', fullfile(pwd,'input','dem','_rawdata'));%, ...
    %'login', login, passwd);


dem.X   = repmat(demTmp.lon, numel(demTmp.lat), 1);
dem.Y   = repmat(demTmp.lat, 1, numel(demTmp.lon));
dem.Y   = flipud(dem.Y);
dem.Z   = double(demTmp.z);
dem.Z   = flipud(dem.Z);
clear demTmp;

% Save
save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');
disp(['DEM saved as: ', 'input/dem/', filename, filesep, filename, '.mat'])

% processSRTM(lat_min, lat_max, lon_min, lon_max, res, filename);

