% Pre-process and run
function runIt(src, ~)

% Retrieve data
part  = guidata(ancestor(src, 'Figure'));
APDTA = getappdata(ancestor(src, 'figure'));  

% Check if input parameters are ok
if part.run_check == 0
    errordlg('There are problems with input parameters, double check.');
    return
end

% If project folder does not exist, then create it
if ~exist(['projects', filesep, part.run_name], 'dir')
    mkdir(['projects', filesep, part.run_name]);
end

% Check if a particle with same name already exists
if exist(['projects', filesep, part.run_name, filesep, part.part.name, '.mat'], 'file')
    choice = questdlg('A particles with the same name already exists. Replace it?', ...
        'Particle', ...
        'Replace', 'Cancel', 'Cancel');
    
    switch choice
        case 'Cancel'
            return
        case 'Replace'
            delete(['projects', filesep, part.run_name, filesep, part.part.name, '.mat']);
            clear_part(src,1,part.part.name);
    end
end

%% RUN 
display('Run started...');
set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', 'Run started, please wait...');
part.traj = get_trajectory(part);
display(part.traj.out_msg); 
set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', '');
display('Done!');

%% Update
update_table(src, part);
APDTA.pltData.(part.part.name) = part;
setappdata(ancestor(src, 'figure'), 'pltData',APDTA.pltData);

% Enable buttons
enableUI(src, 'on');

save(['projects', filesep, part.run_name, filesep, part.part.name, '.mat'], 'part');