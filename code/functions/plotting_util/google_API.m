function google_API(~,~)

if exist(fullfile('code', 'var', 'google_api.mat'),'file') == 2
    tmp = load(fullfile('code', 'var', 'google_api.mat'));
    google_api = tmp.google_api;
else
    google_api = '';
end

apistr = inputdlg('Enter the Google Map API key:', 'API key', [1,50], {google_api});

if isempty(apistr) && isempty(google_api)
    fprintf('The key is empty.\n')
else
    google_api = apistr{1};
    save(fullfile('code', 'var', 'google_api.mat'), 'google_api');
end
