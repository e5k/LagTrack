function dem = download_SRTM(lat_min, lat_max, lon_min, lon_max, res, filename)

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

% Retrieve indices of SRTM files
[lat_minI, lat_maxI, lon_minI, lon_maxI] = get_SRTM_coordinates(lat_min, lat_max, lon_min, lon_max);

% Main 2 loops
for yy = lat_maxI:lat_minI
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
    end   
end

process_SRTM(lat_min, lat_max, lon_min, lon_max, res, filename);

