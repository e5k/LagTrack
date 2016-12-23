% Pre-process and run
function runIt(varargin) %(src, ~)
% Added the varargin option

% Retrieve data
APDTA = getappdata(ancestor(varargin{1}, 'figure'));  

if nargin == 2
    part = {guidata(ancestor(varargin{1}, 'Figure'))};
else
    part = varargin{3};
end

% If project folder does not exist, then create it
if ~exist(['projects', filesep, part{1}.run_name], 'dir')
    mkdir(['projects', filesep, part{1}.run_name]);
end

for iP = 1:length(part)
    % Check if a particle with same name already exists
    if exist(['projects', filesep, part{iP}.run_name, filesep, part{iP}.part.name, '.mat'], 'file')
        choice = questdlg(['A particle named ', part{iP}.part.name, ' already exists in the folder ', part{iP}.run_name, '. Replace it?'], ...
            part{iP}.part.name, ...
            'Replace', 'Cancel', 'Cancel');

        switch choice
            case 'Cancel'
                return
            case 'Replace'
                delete(['projects', filesep, part{iP}.run_name, filesep, part{iP}.part.name, '.mat']);
                clear_part(varargin{1},1,part{iP}.part.name);
        end
    end
end

%% RUN 
set(findobj(ancestor(varargin{1}, 'figure'), 'Tag', 'Errmsg'), 'String', 'Run started, please wait...');
get_trajectory(part);


%% Update
for iP = 1:length(part)
    partTmp     = load(['projects', filesep, part{iP}.run_name, filesep, part{iP}.part.name, '.mat']); % Load particle saved in the get_trajectory function
    update_table(varargin{1}, partTmp.part);                                % Update GUI table
    APDTA.pltData.(partTmp.part.part.name) = partTmp.part;                  % Set appdata
    setappdata(ancestor(varargin{1}, 'figure'), 'pltData',APDTA.pltData);   % Save appdata
    APDTA       = getappdata(ancestor(varargin{1}, 'figure'));              % Retrieve appdata
end
% Enable buttons
set(findobj(ancestor(varargin{1}, 'figure'), 'Tag', 'Errmsg'), 'String', '');
enableUI(varargin{1}, 'on');

%save(['projects', filesep, part.run_name, filesep, part.part.name, '.mat'], 'part');