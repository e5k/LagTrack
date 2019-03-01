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
%           res    :  Resolution (m) (Leave 90 m for no interpolation)
%           name   :  File name, saved in input/dem/
%
%   see also processSRTM, makeDefaultGrid.
% 
% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

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
            dem.res     = res;
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
    dem.res     = res;
    %dem.tiles   = tiles;
    %dem.lat     = LT;
    %dem.lon     = LN;
    dem.type    = 'DEM';
    save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');    
end


demTmp  = readhgt([dem.lat_min, dem.lat_max, dem.lon_min, dem.lon_max], 'interp', 'outdir', fullfile(pwd,'input','dem','_rawdata'), 'srtm1');
dem.X   = repmat(demTmp.lon, numel(demTmp.lat), 1);
dem.Y   = repmat(demTmp.lat, 1, numel(demTmp.lon));
dem.Z   = double(demTmp.z);
clear demTmp;

% Save
save(['input/dem/', filename, filesep, filename, '.mat'], 'dem');

disp('Done!')

% % Main loop
% disp('Accessing SRTM server, please wait...')
% maindir = ['input/dem/', filename];
% for i = 1:length(dem.tiles)
%     outdir  = [maindir, filesep, dem.tiles{i}];    % Tmp directory
%     fprintf('   Downloading SRTM tile %s (%d of %d)... \n', tiles{i}, i, length(dem.tiles))
%     
%     % Check if tile is already downloaded
%     if exist(outdir, 'dir')==7 && exist([outdir, filesep, tiles{i}, '.asc'], 'file')
%         continue
%     end
%     
%     % Else download
%     DL_check = 0;
%     while DL_check == 0
%         DL_check    = 1;    % Check download worked
%         tile_check  = 1;    % Check that the tile exists (i.e. fails in the ocean)
%         try
%             websave([maindir, filesep, tiles{i}, '.zip'], ['http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_v41/SRTM_Data_ArcASCII/', tiles{i}, '.zip']);
%         catch ME
%             if strcmp(ME.identifier, 'MATLAB:webservices:HTTP404StatusCodeError')
%                 tile_check = 0;
%             end
%         end
%         
%         if tile_check == 1 % If the tile exists
%             try
%                 unzip([maindir, filesep, tiles{i}, '.zip'], outdir);
%             catch ME
%                 if strcmp(ME.identifier, 'MATLAB:unzip:invalidZipFile')
%                     DL_check = 0;
%                 end
%             end
%         else % If the file doesn't exist, then create a matrix filled with zeros
%             fprintf('      SRTM tile %s does not exist... Creating an empty tile \n', tiles{i})
%             mkdir([maindir, filesep, tiles{i}]);
%             writeDEM([maindir, filesep, tiles{i}, filesep, tiles{i}, '.asc'],...
%                 dem.lon(i), dem.lat(i,1), zeros(6001,6001), 0.00083333333333333);
%         end
%     end   
% end
% 
%processSRTM(lat_min, lat_max, lon_min, lon_max, res, filename);

