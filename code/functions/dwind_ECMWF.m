
function dwind_ECMWF (lat_min, lat_max, lon_min, lon_max, year_min, year_max, month_min, month_max, filename )

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