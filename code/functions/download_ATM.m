function download_ATM(varargin)%(lat_min, lat_max, lon_min, lon_max, year_min, year_max, month_min, month_max, filename, dataset)
% DOWNLOAD_ATM Download atmospheric data from Reanalysis datasets.
%   DOWNLOAD_ATM(lat_min, lat_max, lon_min, lon_max, year_min, year_max, month_min, month_max, filename, dataset)
%       Download data for the specified spatial and temporal extent with
%       the output name filename. The dataset is either 'Interim' for ECMWF
%       Era-Interim or 'Reanalysis2' for NOAA Reanalysis 2 database.
%
%   See also writeECMWFAPIKey, process_ATM.

% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'Minimum latitude (decimal degree, negative in S hemisphere)', 'Maximum latitude (decimal degree, negative in S hemisphere)', 'Minimum longitude (decimal degree, negative in W hemisphere)', 'Maximum longitude (decimal degree, negative in W hemisphere)', 'Start year (yyyy)', 'End year (yyyy)', 'Start month (mm)', 'End month (mm)', 'Name', 'Dataset (Interim or Reanalysis2)'}, 'Download atmospheric data', 1);
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
elseif nargin == 10   
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

% Make output folder
mkdir(['input/wind/', filename])

%% ERA-INTERIM
if strcmp(dataset, 'Interim')
    
    % Work on input coordinates
    if lon_min < 0; lon_min = 360+lon_min; end
    if lon_max < 0; lon_max = 360+lon_max; end

    txt     = fileread('download_ECMWF_tmp.py');
    txt_new = strrep(txt, 'var_date_start', [num2str(year_min), num2str(month_min, '%02d'), '01']);
    txt_new = strrep(txt_new, 'var_date_end', [num2str(year_max), num2str(month_max, '%02d'), num2str(eomday(year_max, month_max),'%02d')]);
    txt_new = strrep(txt_new, 'var_north', num2str(lat_max));
    txt_new = strrep(txt_new, 'var_south', num2str(lat_min));
    txt_new = strrep(txt_new, 'var_west', num2str(lon_min));
    txt_new = strrep(txt_new, 'var_east', num2str(lon_max));
    txt_new = strrep(txt_new, 'var_out', strrep(['input/wind/', filename, filesep, filename, '.nc'], '\', '/'));
    
    fid = fopen('download_ECMWF.py', 'w');
    fprintf(fid, '%s', txt_new);
    fclose(fid);
    
    !python download_ECMWF.py
    
    delete('download_ECMWF.py');

    
    
    
%% NOAA   
else
    % Work on input coordinates
    if lon_min < 0; lon_min = [num2str(360+lon_min),'E']; end
    if lon_max < 0; lon_max = [num2str(360+lon_max),'E']; end
    if lat_min < 0; lat_min = [num2str(abs(lat_min)), 'S']; else lat_min = [num2str(lat_min), 'N']; end 
    if lat_max < 0; lat_max = [num2str(abs(lat_max)), 'S']; else lat_max = [num2str(lat_max), 'N']; end 
    
    month_list  = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    
    % Reanalysis 1, not used anymore
    if strcmp(dataset, 'Reanalysis1')
        display('Connecting to NOAA NCEP/NCAR Reanalysis 1')
        % Find tid variable
        pge = webread('http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Pressure+Level&Variable=Geopotential+height&group=0&submit=Search');
        fnd = strfind(pge, 'tid=');
        tid = pge(fnd(1)+4:fnd(1)+8);
        
        str1 = 'http://www.esrl.noaa.gov/psd/cgi-bin/DataAccess.pl?DB_dataset=NCEP+Reanalysis+Pressure+Level&DB_variable=Geopotential+height&DB_statistic=Individual+Obs&DB_tid=';
        str2 = '&DB_did=2&DB_vid=14';
        
        % Find nc variable
        pge = webread(strcat(str1, tid, str2));
        fnd = strfind(pge, 'y4.nc');
        nc  = pge(fnd(1)+6:fnd(1)+10);
        
        % Waitbar
        wtb = waitbar(0, '', 'Name', 'Downloading wind...');
        
        %count = 1;
        for i = 1:5
            if i==1
                group = 'gheight';
                name = 'hgt';
                vid = '14';
                unit = 'm';
            elseif i==2
                group = 'uwind';
                name = 'uwnd';
                vid = '18';
                unit = 'm%2Fs';
            elseif i==3
                group = 'vwind';
                name = 'vwnd';
                vid = '19';
                unit = 'm%2Fs';
            elseif i==4
                group = 'relhum';
                name = 'rhum';
                vid = '16';
                unit = '%25';
            elseif i==5
                group = 'temp';
                name = 'air';
                vid = '13';
                unit = 'degK';
            end
            
            for k=1:1*10^9

                page = ['http://www.esrl.noaa.gov/psd/cgi-bin/GrADS.pl?dataset=NCEP+Reanalysis+Pressure+Level&DB_did=2&file=%2FDatasets%2Fncep.reanalysis%2Fpressure%2F',...
                    name,...
                    '.1948.nc+', name,...
                    '.%25y4.nc+', nc,...
                    '&variable=',name,...
                    '&DB_vid=',vid,...
                    '&DB_tid=',tid,...
                    '&units=', unit,...
                    '&longstat=Individual+Obs&DB_statistic=Individual+Obs&stat=&lat-begin=', num2str(lat_min),...
                    '&lat-end=', num2str(lat_max),...
                    '&lon-begin=', num2str(lon_min),...
                    '&lon-end=', num2str(lon_max),...
                    '&dim0=level&level+units=millibar&level=1000.00&level=925.00&level=850.00&level=700.00&level=600.00&level=500.00&level=400.00&level=300.00&level=250.00&level=200.00&level=150.00&level=100.00&level=70.00&level=50.00&level=30.00&level=20.00&level=10.00',...
                    '&dim1=time&year_begin=', num2str(year_min),...
                    '&mon_begin=', month_list{month_min},...
                    '&day_begin=', '1',...
                    '&hour_begin=00+Z',...
                    '&year_end=', num2str(year_max),...
                    '&mon_end=', month_list{month_max},...
                    '&day_end=', num2str(eomday(year_max, month_max)),...
                    '&hour_end=18+Z',...
                    '&X=lon&Y=lat&output=file&bckgrnd=black&use_color=on&cint=&range1=&range2=&scale=100&submit=Create+Plot+or+Subset+of+Data'];
                
                [content, statut] = webread(page);
                if statut == 0
                    break
                else
                    display(['Downloading ', group, '...']) 
                    % Setup  ftp access
                    folder  = 'Public/www/';
                    firs    = strfind(content,'ftp.cdc.noaa.gov/Public/www/');
                    last    = strfind(content,'>FTP a copy of the file');
                    url1    = content(1,firs:last-1);
                    mid     = strfind(url1, 'X');
                    url     = url1(1,mid:length(url1));
                    ftpobj  = ftp('ftp.cdc.noaa.gov');
                    cd(ftpobj, folder);
                    mget(ftpobj, url, '.');
                    
                    % Move files
                    pth     = ['input/wind/', filename, filesep, filename, '_', group, '.nc'];
                    movefile(url, pth);
                    
                    % Updates Waitbar
                    waitbar(i/5, wtb);
                end
                
                if statut == 1
                    break
                end
            end
        end
                
   % Reanalysis 2    
   elseif strcmp(dataset, 'Reanalysis2')
       disp('Connecting to NOAA NCEP/NCAR Reanalysis 2')
        
        % Waitbar
        wtb = waitbar(0, '', 'Name', 'Downloading wind...');
        
        %count = 1;
        for i = 1:5
            if i==1
                group   = 'gheight';
                name    = 'hgt';
                vid     = '1454';
                unit    = 'm';
                var     = 'Geopotential+height';
            elseif i==2
                group   = 'uwind';
                name    = 'uwnd';
                vid     = '4300';
                unit    = 'm%2Fs';
                var     = 'U-wind';
            elseif i==3
                group   = 'vwind';
                name    = 'vwnd';
                vid     = '4306';
                unit    = 'm%2Fs';
                var     = 'V-wind';
            elseif i==4
                group   = 'relhum';
                name    = 'rhum';
                vid     = '4271';
                unit    = '%25';
                var     = 'Relative+Humidity';
            elseif i==5
                group   = 'temp';
                name    = 'air';
                vid     = '1447';
                unit    = 'degK';
                var     = 'Air+Temperature';
            end
            
            [tid, nc] = get_tid(var, vid);
            
          %  for k=1:1*10^9
                
                page = ['http://www.esrl.noaa.gov/psd/cgi-bin/GrADS.pl?dataset=NCEP%2FDOE+AMIP-II+Reanalysis+%28Reanalysis-2%29&DB_did=61&file=%2FDatasets%2Fncep.reanalysis2%2Fpressure%2F',...
                    name,...
                    '.1979.nc+', name,...
                    '.%25y4.nc+', nc,...
                    '&variable=',name,...
                    '&DB_vid=',vid,...
                    '&DB_tid=',tid,...
                    '&units=', unit,...
                    '&longstat=Individual+Obs&DB_statistic=Individual+Obs&stat=&lat-begin=', num2str(lat_min),...
                    '&lat-end=', num2str(lat_max),...
                    '&lon-begin=', num2str(lon_min),...
                    '&lon-end=', num2str(lon_max),...
                    '&dim0=level&level+units=millibar&level=1000.00&level=925.00&level=850.00&level=700.00&level=600.00&level=500.00&level=400.00&level=300.00&level=250.00&level=200.00&level=150.00&level=100.00&level=70.00&level=50.00&level=30.00&level=20.00&level=10.00',...
                    '&dim1=time&year_begin=', num2str(year_min),...
                    '&mon_begin=', month_list{month_min},...
                    '&day_begin=', '1',...
                    '&hour_begin=00+Z',...
                    '&year_end=', num2str(year_max),...
                    '&mon_end=', month_list{month_max},...
                    '&day_end=', num2str(eomday(year_max, month_max)),...
                    '&hour_end=18+Z',...
                    '&X=lon&Y=lat&output=file&bckgrnd=black&use_color=on&fill=lines&cint=&range1=&range2=&scale=100&maskf=%2FDatasets%2Fncep.reanalysis2%2Fgaussian_grid%2Fland.sfc.gauss.nc&maskv=Land-sea+mask&submit=Create+Plot+or+Subset+of+Data'];
          
%                 [content, statut] = urlread(page);
%                 if statut == 0
%                     break
%                 else

                    content = webread(page, 'timeout', inf);
                    disp(['Downloading ', group, '...'])
                    % Setup  ftp access
                    folder  = 'Public/www/';
                    firs    = strfind(content,'ftp.cdc.noaa.gov/Public/www/');
                    last    = strfind(content,'>FTP a copy of the file');
                    url1    = content(1,firs:last-1);
                    mid     = strfind(url1, 'X');
                    url     = url1(1,mid:length(url1));
                    ftpobj  = ftp('ftp.cdc.noaa.gov');
                    cd(ftpobj, folder);
                    mget(ftpobj, url, '.');
                    
                    % Move files
                    pth     = ['input/wind/', filename, filesep, filename, '_', group, '.nc'];
                    movefile(url, pth);
                    
                    % Updates Waitbar
                    waitbar(i/5, wtb);
%                end
                
%                 if statut == 1
%                     break
%                 end
%            end
        end       
    end    
end

process_ATM(filename, dataset)


function [tid, nc] = get_tid(var, vid)

%pge = urlread(['http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP/DOE+AMIP-II+Reanalysis+(Reanalysis-2)+&Variable=', var]);
pge = webread(['http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?&Dataset=NCEP/DOE+AMIP-II+Reanalysis+(Reanalysis-2)+&Variable=', var]);

fnd = strfind(pge, 'tid=');
tid = pge(fnd(1)+4:fnd(1)+8);

pge2 = webread(['http://www.esrl.noaa.gov/psd/cgi-bin/DataAccess.pl?DB_dataset=NCEP/DOE+AMIP-II+Reanalysis+(Reanalysis-2)&DB_variable=', var,...
    '&DB_statistic=Individual+Obs&DB_tid=', tid,...
    '&DB_did=61&DB_vid=', num2str(vid)]);
fnd = strfind(pge2, 'y4.nc ');
nc  = pge2(fnd(1)+6:fnd(1)+10);
        