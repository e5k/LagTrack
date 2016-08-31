% Clear particles

function delete_part(src,~)

choice = questdlg('Do you want to permanently delete these particles?', ...
	'Delete particles', ...
	'Yes','No','No');
% Handle response
switch choice
    case 'Yes'
    case 'No'; return
end

% Load GUI data
APDTA   = getappdata(ancestor(src, 'figure'));   
fld     = fieldnames(APDTA.pltData);
pltData = APDTA.pltData;

% Get table data from GUI
List    = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String');
Tab     = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataTable'), 'Data');
ListV   = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'Value');

% Delete files
for i = 1:length(ListV)
    tmp = pltData.(List{ListV(i)});
    fl  = ['projects', filesep, tmp.run_name, filesep, tmp.part.name, '.mat'];
    if exist(fl, 'file')
        delete(fl);
    else
        error(['File ', fl, ' not found']);
    end
end

% Setup logical vector
partCheck = zeros(length(List),1);
if exist('ListV', 'var')
    partCheck(ListV) = true;
end

% Remove selected data from GUI data
for i = 1:length(fld)
    if partCheck(i) == 1
        APDTA.pltData = rmfield(APDTA.pltData, fld{i});
    end
end

% Remove entries
List(logical(partCheck)) = [];
Tab(logical(partCheck),:) = [];

% Update table data
set(findobj(ancestor(src, 'figure'), 'Tag', 'DataTable'), 'Data', Tab);
set(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'Value', []);
set(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String', List);

% If all particles were cleared, disable interface
if isempty(List)
    enableUI(src, 'off');
end

% Update GUI data
setappdata(ancestor(src, 'figure'), 'pltData', APDTA.pltData);   