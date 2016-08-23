% Clear particles

function clear_part(src,~)

% Load GUI data
APDTA = getappdata(ancestor(src, 'figure'));   
fld   = fieldnames(APDTA.pltData)

% Get table data from GUI
List    = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String');
ListV   = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'Value');
Tab     = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataTable'), 'Data');

% Setup logical vector
partCheck = zeros(length(List),1);
partCheck(ListV) = true;

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