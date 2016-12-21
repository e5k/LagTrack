function makeDefaultGrid(varargin)

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'Mean altitude (m asl)', 'Name'}, ' ', 1);
    if isempty(answer)
        return
    end
    alt         = str2double(answer{2});
    grdName     = answer{2};    
elseif nargin == 3   
    alt         = varargin{1};
    grdName     = varargin{3};
else
    error('Wrong number of input arguments, should be either 0 or 3');
end

grdName = [grdName, '_STD'];

% Check if folder already exist
if exist(['input/dem/', grdName], 'dir') == 7
    choice = questdlg('A folder with the same name already exists. Do you want to delete it?', ...
        'Wind Name', ...
        'No', 'Yes', 'Yes');
    % Handle response
    switch choice
        case 'Yes'
            rmdir(['input/dem/', grdName], 's');
        case 'No'
            return
    end
end

% Make output folder
mkdir(['input/dem/', grdName])

dem.lat_min = 0;
dem.lat_max = 0;
dem.lon_min = 0;
dem.lon_max = 0;
dem.res     = 0;
dem.type    = 'GRID';
dem.X       = 0;
dem.Y       = 0;
dem.Z       = alt;

% Save data
save(['input/dem/', grdName, filesep, grdName, '.mat'], 'dem');
disp(['Saved as ', 'input/wind/', grdName, filesep, grdName, '.mat']);