function writeAPI(~,~,type)
% Type: 0 = ERA-Interim
%       1 = ERA-5
% writeECMWFAPIKey Writes ECMWF API key to home folder.
%
%   See also downloadATM, preprocessATM.

% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

if ispc
    userdir= getenv('USERPROFILE'); 
else
    userdir= getenv('HOME');
end

% Define if ERA-Interim or ERA-5
if type == 0
    targetFile = '.ecmwfapirc';
    defStr = {sprintf('{\n\t"url"   : "https://api.ecmwf.int/v1",\n\t"key"   : "___your ID___",\n\t"email" : "___your email___"\n}')};
else
    targetFile = '.cdsapirc';
    defStr = {sprintf('url: https://cds.climate.copernicus.eu/api/v2\nkey: ___yourUID___:___yourAPIKey___')};
end

% Check if exists
if exist([userdir, filesep, targetFile], 'file')
    choice = questdlg('It seems that an API key already exists. Overwrite?', ...
    'API key', ...
    'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            choice = 1;
        case 'No'
            choice = 0;
    end
else
    choice = 1;
end

% Write the key
if choice == 1
    apistr = inputdlg('Enter the content of the API key:', 'API key', [5,100], defStr);
    if isempty(apistr)
        return
    end
    
    apistr = apistr{1};

    fid = fopen([userdir, filesep, targetFile], 'w');
    for i = 1:size(apistr,1)
        for j = 1:size(apistr,2)

            if j == size(apistr,2)
                fprintf(fid, '%s\n', apistr(i,j));
            else
                fprintf(fid, '%s', apistr(i,j));
            end
        end
    end
    fclose(fid);
end