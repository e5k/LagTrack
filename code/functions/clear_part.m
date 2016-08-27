% Clear particles

function clear_part(src,~,varargin)

% Load GUI data
APDTA = getappdata(ancestor(src, 'figure'));   
fld   = fieldnames(APDTA.pltData);

% Get table data from GUI
List    = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String');
Tab     = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataTable'), 'Data');

if nargin == 2    % If this function is called from the main GUI
    % Get table data from GUI
    ListV   = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'Value');

else % Function called when starting the run
    for i = 1:length(List)
        if strcmp(varargin{1}, List{i}); ListV = i; end           
    end
end
    
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