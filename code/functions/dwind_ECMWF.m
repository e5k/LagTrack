
function dwind_ECMWF (lat_min, lat_max, lon_min, lon_max, year_min, year_max, month_min, month_max, filename, dataset )

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
    
    if lat_min < 0; lat_min = 360+lat_min; end
    if lat_max < 0; lat_min = 360+lat_max; end
    
    month_list  = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    
    
    if strcmp(dataset, 'Reanalysis1')
        % Find tid variable
        pge = urlread('http://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis+Pressure+Level&Variable=Geopotential+height&group=0&submit=Search');
        fnd = strfind(pge, 'tid=');
        tid = pge(fnd(1)+4:fnd(1)+8);
        
        str1 = 'http://www.esrl.noaa.gov/psd/cgi-bin/DataAccess.pl?DB_dataset=NCEP+Reanalysis+Pressure+Level&DB_variable=Geopotential+height&DB_statistic=Individual+Obs&DB_tid=';
        str2 = '&DB_did=2&DB_vid=14';
        
        % Find nc variable
        pge = urlread(strcat(str1, tid, str2));
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
                
                %                             http://www.esrl.noaa.gov/psd/cgi-bin/GrADS.pl?dataset=NCEP%2FDOE+AMIP-II+Reanalysis+%28Reanalysis-2%29&DB_did=61&file=%2FDatasets%2Fncep.reanalysis2%2Fpressure%2F...
                %                             .1979.nc+
                %                             .%25y4.nc+
                %                             &variable=
                %                             &DB_vid=
                %                             &DB_tid=
                %                             &units=
                %                             &longstat=Individual+Obs&DB_statistic=Individual+Obs&stat=&lat-begin=
                %                             &lat-end=
                %                             &lon-begin=
                %                             &lon-end=
                %                             &dim0=level&level+units=millibar&level=1000.00&level=925.00&level=850.00&level=700.00&level=600.00&level=500.00&level=400.00&level=300.00&level=250.00&level=200.00&level=150.00&level=100.00&level=70.00&level=50.00&level=30.00&level=20.00&level=10.00
                %                             &dim1=time&year_begin=
                %                             &mon_begin=
                %                             &day_begin=
                %                             &hour_begin=
                %                             +Z&year_end=
                %                             &mon_end=
                %                             &day_end=
                %                             &hour_end=
                %                             +Z&X=lon&Y=lat&output=file&bckgrnd=black&use_color=on&fill=lines&cint=&range1=&range2=&scale=100&maskf=%2FDatasets%2Fncep.reanalysis2%2Fgaussian_grid%2Fland.sfc.gauss.nc&maskv=Land-sea+mask&submit=Create+Plot+or+Subset+of+Data
                
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
                
                [content, statut] = urlread(page);
                if statut == 0
                    break
                else
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
                    delete(url)
                    
                    % Updates Waitbar
                    waitbar(count/5, wtb);
                    count = count+1;
                end
                
                if statut == 1;
                    break
                end
            end
        end
        
        
        
        
        
        
    elseif strcmp(dataset, 'Reanalysis2')
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        for i = group_start:3
            k=1;
            if i==1
                group = 'gheight';
                name = 'hgt';
                vid = '1454';
                unit = 'm';
            elseif i==2
                group = 'uwind';
                name = 'uwnd';
                vid = '4300';
                unit = 'm%2Fs';
            elseif i==3
                group = 'vwind';
                name = 'vwnd';
                vid = '4306';
                unit = 'm%2Fs';
            elseif i==4
                group = 'relhum';
                name = 'rhum';
                vid = '4271';
                unit = '%25';
            elseif i==5
                group = 'temp';
                name = 'air';
                vid = '1447';
                unit = 'degK';
            end
            
        end
        
        
    end
    
end

preprocess_ECMWF(filename)