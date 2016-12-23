function atmName = makeStandardAtm(varargin)
% makeStandardAtm Create an atmosphere file based on the 1976 Standard
% atmosphere specified in code/var/1976_standard_ATM.xlsx with a custom
% wind velocity constant with height.
%   makeStandardAtm
%       Open the GUI
%   makeStandardAtm(uwind, vwind, atmName)
%       Download the SRTM tiles covering the specified extent
%           uwind   :  U wind velocity (m/s)
%           vwind   :  V wind velocity (m/s)
%           atmName :  File name, saved in input/wind/
%
%   see also downloadSRTM, makeStandardAtm.
% 
% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'U wind (m/s, positive towards E)', 'V wind (m/s, positive towards N)', 'Name'}, ' ', 1);
    if isempty(answer)
        return
    end
    uwind       = str2double(answer{1});
    vwind       = str2double(answer{2});
    atmName     = answer{3};    
elseif nargin == 3   
    uwind       = varargin{1};
    vwind       = varargin{2};
    atmName     = varargin{3};
else
    error('Wrong number of input arguments, should be either 0 or 3');
end

atmName = [atmName, '_STD'];

% Check if folder already exist
if exist(['input/wind/', atmName], 'dir') == 7
    choice = questdlg('A folder with the same name already exists. Do you want to delete it?', ...
        'Wind Name', ...
        'No', 'Yes', 'Yes');
    % Handle response
    switch choice
        case 'Yes'
            rmdir(['input/wind/', atmName], 's');
        case 'No'
            return
    end
end

% Make output folder
mkdir(['input/wind/', atmName])



% Retrieve the standard atmosphere
a = xlsread('1976_standard_ATM.xlsx');
rhoair = zeros(1,1, size(a,1), 2);
rhoair(:,:,:,1) = reshape(a(:,2), 1,1,size(a,1),1);
rhoair(:,:,:,2) = reshape(a(:,2), 1,1,size(a,1),1);

muair = zeros(1,1, size(a,1), 2);
muair(:,:,:,1) = reshape(a(:,3), 1,1,size(a,1),1);
muair(:,:,:,2) = reshape(a(:,3), 1,1,size(a,1),1);

alt = zeros(1,1, size(a,1), 2);
alt(:,:,:,1) = reshape(a(:,1), 1,1,size(a,1),1);
alt(:,:,:,2) = reshape(a(:,1), 1,1,size(a,1),1);

u = zeros(1,1, size(a,1), 2);
u(:,:,:,1) = ones(1,1,size(a,1),1).*uwind;
u(:,:,:,2) = ones(1,1,size(a,1),1).*uwind;

v = zeros(1,1, size(a,1), 2);
v(:,:,:,1) = ones(1,1,size(a,1),1).*vwind;
v(:,:,:,2) = ones(1,1,size(a,1),1).*vwind;

time = [datenum([1948,1,1]); now];

atm.lat = 0;
atm.lon = 0;
atm.time = time;
atm.alt = alt;
atm.u = u;
atm.v = v;
atm.rhoair = rhoair;
atm.muair = muair;

% Save data
save(['input/wind/', atmName, filesep, atmName, '.mat'], 'atm');
disp(['Saved as ', 'input/wind/', atmName, filesep, atmName, '.mat']);