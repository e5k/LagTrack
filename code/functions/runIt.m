% Pre-process and run
function runIt(src, ~)

% Retrieve data
part = guidata(ancestor(src, 'Figure'));

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
            rm(['projects', filesep, part.run_name, filesep, part.part.name, '.mat']);
    end
end

