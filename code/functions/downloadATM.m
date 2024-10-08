function downloadATM(varargin)
% downloadATM Download atmospheric data from Reanalysis datasets.
%   downloadATM
%       Opens the GUI to download atmospheric data from Reanalysis datasets
%   downloadATM(lat_min, lat_max, lon_min, lon_max, year_min, year_max, month_min, month_max, filename, dataset)
%       Download data for the specified spatial and temporal extent with
%       the output name filename. The dataset is either 'Interim' for ECMWF
%       Era-Interim, 'ERA5', 'Reanalysis1' for NOAA Reanalysis 1 and 'Reanalysis2' 
%       for NOAA Reanalysis 2. When using ERA5, it is also possible to add
%       the number of hours between wind profiles as a last argument:
%   downloadATM(lat_min, lat_max, lon_min, lon_max, year_min, year_max, month_min, month_max, filename, 'ERA5', nHours)
%       nHours accepts 1, 2, 3, 4, 6, 8 and 12. Default is 6.
%
%   See also processATM, displayATM, makeStandardAtm.
%
% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

% variables initiation
nHours = 6; % Number of hours between wind profiles. Only used for ERA-5

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'Minimum latitude (decimal degree, negative in S hemisphere)', 'Maximum latitude (decimal degree, negative in S hemisphere)', 'Minimum longitude (decimal degree, negative in W hemisphere)', 'Maximum longitude (decimal degree, negative in W hemisphere)', 'Start year (yyyy)', 'End year (yyyy)', 'Start month (mm)', 'End month (mm)', 'Name', 'Dataset (ERA5, Interim, Reanalysis1 or Reanalysis2)'}, 'Download atmospheric data', 1);
    if isempty(answer)
        return
    end
    lat_min     = str2double(answer{1});
    lat_max     = str2double(answer{2});
    lon_min     = str2double(answer{3});   
    lon_max     = str2double(answer{4});
    year_min    = str2double(answer{5});
    year_max    = str2double(answer{6});
    month_min   = str2double(answer{7});
    month_max   = str2double(answer{8});
    filename    = answer{9};
    dataset     = answer{10};
elseif nargin == 10 || nargin == 11
    lat_min     = varargin{1};
    lat_max     = varargin{2};
    lon_min     = varargin{3};   
    lon_max     = varargin{4};
    year_min    = varargin{5};
    year_max    = varargin{6};
    month_min   = varargin{7};
    month_max   = varargin{8};
    filename    = varargin{9};
    dataset     = varargin{10};
    % If ERA 5 and last argument
    if strcmp(dataset, 'ERA5') && nargin == 11
        nHours = varargin{11};
        if ~ismember(nHours, [1,2,3,4,6,8,12,24])
            error('nHours accepts 1, 2, 3, 4, 6, 8, 12 or 24')
        end        
    end
else
    error('Wrong number of input arguments, should be either 0 or 10');
end


% Extend the range
if abs(lat_max - lat_min) < 2.5; lat_min = lat_min-1.5; lat_max = lat_max+1.5; end
if abs(lon_max - lon_min) < 2.5; lon_min = lon_min-1.5; lon_max = lon_max+1.5; end

% Check if folder already exist
if exist(['input/wind/', filename], 'dir') == 7
    choice = questdlg('A folder with the same name already exists. Do you want to delete it?', ...
        'Wind Name', ...
        'No', 'Yes', 'Yes');
    % Handle response
    switch choice
        case 'Yes'
            rmdir(['input/wind/', filename], 's');
        case 'No'
            return
    end
end

% If ERA-5, prepare the hours vector
if strcmp(dataset, 'ERA5')
    vHours = 0:nHours:23;
    vHours = datestr(datetime([ones(size(vHours))', ones(size(vHours))', ones(size(vHours))', vHours', zeros(size(vHours))', zeros(size(vHours))']), 15);
    t = '[';
    for i=1:size(vHours,1)
        t = [t, '''',vHours(i,:),''', '];
    end
    t = [t(1:end-2), ']'];
end

% Make output folder
mkdir(['input/wind/', filename])

%% ECMWF
if strcmp(dataset, 'Interim') || strcmp(dataset, 'ERA5') 
    
    % Work on input coordinates
    if lon_min < 0; lon_min = 360+lon_min; end
    if lon_max < 0; lon_max = 360+lon_max; end
    
    if strcmp(dataset, 'Interim')
        fl2read = 'code/functions/dependencies/ecmwf-api-client-python/download_Interim_tmp.py';
    else
        fl2read = 'code/functions/dependencies/cdsapi-0.2.5/download_ERA5_tmp.py';
    end
    
    txt     = fileread(fl2read);
    txt_new = strrep(txt, 'var_year_start', num2str(year_min));
    txt_new = strrep(txt_new, 'var_year_end', num2str(year_max));
    txt_new = strrep(txt_new, 'var_month_start', num2str(month_min));
    txt_new = strrep(txt_new, 'var_month_end', num2str(month_max));
    txt_new = strrep(txt_new, 'var_north', num2str(lat_max));
    txt_new = strrep(txt_new, 'var_south', num2str(lat_min));
    txt_new = strrep(txt_new, 'var_west', num2str(lon_min));
    txt_new = strrep(txt_new, 'var_east', num2str(lon_max));
    txt_new = strrep(txt_new, 'var_out', strrep(['input/wind/', filename, filesep, filename], '\', '/'));
    
    if strcmp(dataset, 'ERA5')
        txt_new = strrep(txt_new, 'nHours', t);
    end
    
    fid = fopen('download_ECMWF.py', 'w');
    fprintf(fid, '%s', txt_new);
    fclose(fid);
    
    !python download_ECMWF.py
    
    %delete('download_ECMWF.py');

%% NOAA
else
    % Work on input coordinates
    if lon_min < 0; lon_min = 360+lon_min; end
    if lon_max < 0; lon_max = 360+lon_max; end
    
    % Reanalysis
    %  - Now downloads worldwide files per year, and the zone extration is
    %    done during post processing. The download time is longer, but files
    %    are preserved, so donwload is skipped if the file already exists.
    %  - Now works the same for Reanalysis 1 and 2
    
    
    if strcmp(dataset, 'Reanalysis1')
        target_dir = 'input/wind/_Reanalysis1_Rawdata/';
        ftp_dir    = 'Datasets/ncep.reanalysis/pressure/';
    elseif strcmp(dataset, 'Reanalysis2')
        target_dir = 'input/wind/_Reanalysis2_Rawdata/';
        ftp_dir    = 'Datasets/ncep.reanalysis2/pressure/';
    else
        error('Unknown dataset requested')
    end
    
    disp('Connecting to NOAA...')
    
    % Files are now preserved in the folder Reanalyisi_data/
    if ~exist(target_dir, 'dir')
        mkdir(target_dir);
    end
    
    varList = {'hgt', 'uwnd', 'vwnd', 'rhum', 'air'}; % List of variables to download
    
    for iV = 1:length(varList)
        for iY = year_min:year_max
            fl = [varList{iV}, '.', num2str(iY), '.nc'];
            if exist([target_dir, fl], 'file')     % If the file exists
                fprintf('\t%s already exists, skipping download...\n', fl)
            else                                                        % Else request ftp
                fprintf('\tDownloading %s, please wait...\n', fl)
                ftpobj  = ftp('ftp.cdc.noaa.gov');
                cd(ftpobj, ftp_dir);
                mget(ftpobj, fl, target_dir);
            end
        end
        
    end
end
processATM(filename, dataset, lat_min, lat_max, lon_min, lon_max, year_min, year_max, month_min, month_max)
