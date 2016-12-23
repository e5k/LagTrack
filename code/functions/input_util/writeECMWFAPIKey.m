function writeECMWFAPIKey
% writeECMWFAPIKey Writes ECMWF API key to home folder.
%
%   See also downloadATM, preprocessATM.

% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

if ispc; 
    userdir= getenv('USERPROFILE'); 
else
    userdir= getenv('HOME');
end

if exist([userdir, filesep, '.ecmwfapirc'], 'file')
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

if choice == 1
    apistr = inputdlg('Enter the content of the API key:', 'API key', [5,100]);
    apistr = apistr{1};

    fid = fopen([userdir, filesep, '.ecmwfapirc'], 'w');
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